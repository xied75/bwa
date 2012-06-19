#include <Windows.h>
#include <time.h>
#include <malloc.h>
#include <stdio.h>
#include "bwtaln.h"
#include "bwase.h"
#include "bntseq.h"
#include "kvec.h"
#include "khash.h"
#include "utils.h"
#include "glibc_win64_flat\stdlib.h"

#include "win64_util\Win32_util.h"

#ifdef THREAD

#define READ_BATCH 0x40000

//MT
DWORD WINAPI pe_worker(LPVOID);
HANDLE *pe_eid_a; //HANDLE array for event objects, type a, to serialize access to fq files
HANDLE *pe_eid_b; //HANDLE array for event objects, type b, to serialize fprintf to SAM file

typedef struct {
    int t;
    int n_threads;
    int batch;
    LPVOID lp_sa[2]; //to hold the mmf pointers
    int * sai_i[2]; //to hold the index of SAIs
    bwa_seqio_t *ks[2];
    gap_opt_t opt0;
    bntseq_t *bns;
    const char *prefix;
    bwt_t *bwt;
    pe_opt_t *popt;
    uint8_t *pac;
    bntseq_t *ntbns;
} pe_thread_data_t; //to pass data into threads


int tot_seqs = 0; //to be visible to all threads thus moved here, originally in bwa_sai2sam_pe_core()
//END MT


typedef struct {
    int n;
    bwtint_t *a;
} poslist_t;

typedef struct {
    uint64_t x, y;
} b128_t;

#define b128_lt(a, b) ((a).x < (b).x)
#define b128_eq(a, b) ((a).x == (b).x && (a).y == (b).y)
#define b128_hash(a) ((uint32_t)(a).x)

KHASH_INIT(b128, b128_t, poslist_t, 1, b128_hash, b128_eq)

typedef struct {
    kvec_t(b128_t) arr;
    kvec_t(b128_t) pos[2];
    kvec_t(bwt_aln1_t) aln[2];
} pe_data_t;

#define MIN_HASH_WIDTH 1000

typedef struct {
    double avg, std, ap_prior;
    bwtint_t low, high, high_bayesian;
} isize_info_t;


extern int g_log_n[256]; // in bwase.c
static kh_b128_t *g_hash;

//in bwape.c
extern bntseq_t *bwa_open_nt(const char *prefix);
extern void bwa_print_sam_SQ(const bntseq_t *bns);
extern void bwa_print_sam_PG();
extern ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii);
extern void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
extern int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, double ap_prior, int64_t L);
extern int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii);

//in bwase.c
extern int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);

//mt code
int * scan_sai(LPVOID lpView);
int bwa_cal_pac_pos_pe_mt(const bntseq_t *bns, const char *prefix, bwt_t *const _bwt, int n_seqs, bwa_seq_t *seqs[2], LPVOID *lp_sa[2], isize_info_t *ii,
                       const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii);
bwt_t *bwt_restore_bwt_mmf(char *fn);
void bwt_restore_sa_mmf(char *fn, bwt_t *bwt);

void bwa_sai2sam_pe_core_mt(const char *prefix, char *const fn_sa[2], char *const fn_fa[2], pe_opt_t *popt)
{
    extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
    int i, j = 0;
    bwa_seqio_t *ks[2];
    bntseq_t *bns, *ntbns = 0;

    LPVOID lp_sa[2]; //to hold the pointer to memory mapped files: SAI1, SAI2
    int *sa_index[2]; // to hold the sai scan results
    int batch = 0; //how many batches to process
    LPVOID lpView_pac; // to hold mmf of PAC

    gap_opt_t opt, opt0;
    khint_t iter;
    char str[1024];
    bwt_t *bwt;
    uint8_t *pac;

    // initialization
    bwase_initialize(); // initialize g_log_n[] in bwase.c
    pac = 0; bwt = 0;
    
    bns = bns_restore(prefix);
    srand48(bns->seed);
    lp_sa[0] = mmf_open(fn_sa[0]);
    lp_sa[1] = mmf_open(fn_sa[1]);

    g_hash = kh_init(b128);
    
    memcpy(&opt, (uint32_t*)lp_sa[0], sizeof(gap_opt_t));

    ks[0] = bwa_open_reads(opt.mode, fn_fa[0]);
    opt0 = opt;
    memcpy(&opt, (uint32_t*)lp_sa[1], sizeof(gap_opt_t)); // overwritten!

    ks[1] = bwa_open_reads(opt.mode, fn_fa[1]);
    if (!(opt.mode & BWA_MODE_COMPREAD)) {
        popt->type = BWA_PET_SOLID;
        ntbns = bwa_open_nt(prefix);
    } else { // for Illumina alignment only
        if (popt->is_preload) {
            HANDLE hProc; BOOL b;
            hProc = GetCurrentProcess();
            b = SetProcessWorkingSetSize(hProc, 5000000000, 6000000000);
            if (!b) fprintf(stderr, "Failed to set new size\n");

            strcpy(str, prefix); strcat(str, ".pac");
            lpView_pac = mmf_open(str);
            pac = (uint8_t*)lpView_pac;

            strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt_mmf(str);
            strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa_mmf(str, bwt);
            /*pac = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
            rewind(bns->fp_pac);
            fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);*/
        }
    }

    sa_index[0] = scan_sai(lp_sa[0]);
    sa_index[1] = scan_sai(lp_sa[1]);
    batch = _msize(sa_index[0])/sizeof(int);

    bwa_print_sam_SQ(bns);
    bwa_print_sam_PG();

    if (popt->n_threads > batch) popt->n_threads = batch; //if you have more cores than batch, use batch
    
    {//Multithreading block
        DWORD ThreadID;
        HANDLE *tid; //HANDLE array for thread objects

        int j = 0;

        pe_thread_data_t * pe_thread_data = (pe_thread_data_t*)calloc(popt->n_threads, sizeof(pe_thread_data_t));
        pe_eid_a = (HANDLE*)calloc(popt->n_threads, sizeof(HANDLE));
        pe_eid_b = (HANDLE*)calloc(popt->n_threads, sizeof(HANDLE));
        tid = (HANDLE*)calloc(popt->n_threads, sizeof(HANDLE));

        //create event objects
        for(; j < popt->n_threads; j++) {
            pe_eid_a[j] = CreateEvent(NULL, FALSE, FALSE, NULL);
            if (pe_eid_a[j] == NULL) 
            { 
                printf("CreateEvent failed (%d)\n", GetLastError()); //To Change, don't leave like this!!!
                return;
            }
            pe_eid_b[j] = CreateEvent(NULL, FALSE, FALSE, NULL);
            if (pe_eid_b[j] == NULL) 
            { 
                printf("CreateEvent failed (%d)\n", GetLastError()); //To Change, don't leave like this!!!
                return;
            }
        }

        //create threads
        for(j=0; j < popt->n_threads; j++) {
            pe_thread_data[j].t = j;
            pe_thread_data[j].batch = batch;
            pe_thread_data[j].n_threads = popt->n_threads;
            pe_thread_data[j].lp_sa[0] = lp_sa[0];
            pe_thread_data[j].lp_sa[1] = lp_sa[1];
            pe_thread_data[j].sai_i[0] = sa_index[0];
            pe_thread_data[j].sai_i[1] = sa_index[1];
            pe_thread_data[j].ks[0] = ks[0];
            pe_thread_data[j].ks[1] = ks[1];
            pe_thread_data[j].opt0 = opt0;
            pe_thread_data[j].bns = bns;
            pe_thread_data[j].prefix = prefix;
            pe_thread_data[j].bwt = bwt;
            pe_thread_data[j].popt = popt;
            pe_thread_data[j].pac = pac;
            pe_thread_data[j].ntbns = ntbns;

            tid[j] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) pe_worker, &pe_thread_data[j], 0, &ThreadID);
            if (tid[j] == NULL)
            {
                printf("CreateThread error: %d\n", GetLastError()); // To change, don't leave like this!!!
                return;
            }
        }

        WaitForMultipleObjects(popt->n_threads, &tid[0], TRUE, INFINITE);

        free(pe_thread_data);
        free(pe_eid_a);
        free(pe_eid_b);
        free(tid);
    }//END of Multithreading block

    // destroy
    bns_destroy(bns);
    if (ntbns) bns_destroy(ntbns);
    for (i = 0; i < 2; ++i) {
        bwa_seq_close(ks[i]);
        
        //fclose(fp_sa[i]);
        mmf_close(lp_sa[i]);
        free(sa_index[i]);
    }
    for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
        if (kh_exist(g_hash, iter)) free(kh_val(g_hash, iter).a);
    kh_destroy(b128, g_hash);
    if (pac) {
        free(pac); bwt_destroy(bwt); //to change, call mmf_close(), with adjusted address maybe
        // on two mmfs
        mmf_close(lpView_pac); //free PAC mmf
    }
}

int * scan_sai(LPVOID lpView)
{
    MEMORY_BASIC_INFORMATION MBI;
    int n_aln, offset = 0, count = 0, max, i = 0;
    int *keep; //assume SAI file length divide by 4 won't exceed int range!

    const int keep_init = 5000; //will take 5000*4=20k mem, assume the biggest SAI/fq has 1 billion reads
                                //divide by 0x40000 is < 5000, doing this to avoid checking alloc inside loop
    const int header = 16; // 64/4, sai header 64bytes, we cast address to uint32 thus 4
    const int t_step = sizeof(bwt_aln1_t)/4;
    keep = (int*)calloc(keep_init, sizeof(int)); 

    VirtualQuery((LPCVOID)lpView, &MBI, sizeof(MBI));
    max = MBI.RegionSize/4;

    offset = header; //the first n_aln;

    while(offset < max)
    {
        n_aln = *((uint32_t*)lpView+offset);
        if (count % READ_BATCH == 0) keep[i++] = offset; //since count starts from 0, here we ARE taking the NEXT offset
                                                         //also 0 % READ_BATCH == 0, thus keep[0] always=16
        if (n_aln == 0)
            offset++;
        else {
            offset += n_aln*t_step+1;
        }
        count++;
    }

    //trim keep to the real size, no need to add 1 on i
    keep = (int*)realloc(keep, i*sizeof(int));

    return keep; //please remember to free() this! Beware that MBI.RegionSize is rounded up to 4kB page size step thus
                 // >= real file size, so when using keep, keep[x] might not be real range since file ended after keep[x-1].
}

DWORD WINAPI pe_worker( LPVOID data ) {

    //threading related
    pe_thread_data_t * pe_data_mt = (pe_thread_data_t*)data;
    int i = 0;
    const int t = pe_data_mt->t; //t is the thread no. e.g. t = 0, (suppose total thread number is 4)
    int b = t; //b will change during loop. e.g. b will be 0, 4, 8, 12, etc.
    int signal_t, loop_count = 1;
    //threading related end

    //BWA related
    int n_seqs = 0;
    bwa_seq_t *seqs[2];
    clock_t ct; //renamed from t because I used t above
    int cnt_chg;
    ubyte_t *pacseq;
    isize_info_t last_ii;
    isize_info_t ii;
    last_ii.avg = -1.0;

    //loop over until all batches for this thread is done
    while (1)
    {
        fprintf(stderr, "T[%u] here, to process batch [%u]\n", t, b);
        
        fprintf(stderr, "T[%u] going to wait on event_a\n", t);

        //wait for event a, except t0&b0, i.e. the first batch
        if (t != 0 || b != 0)
            WaitForSingleObject(pe_eid_a[t], INFINITE);

        //works1, serialized area
        /*fprintf(stderr, "T[%u] working_1 on batch [%u]\n", t, b);
        for(i=0;i<2;i++){
            fprintf(stderr, "[%u+%u]", t, b);
        }*/

        seqs[0] = bwa_read_seq(pe_data_mt->ks[0], 0x40000, &n_seqs, pe_data_mt->opt0.mode, pe_data_mt->opt0.trim_qual);
        if (seqs[0] != 0) {
            seqs[1] = bwa_read_seq(pe_data_mt->ks[1], 0x40000, &n_seqs, pe_data_mt->opt0.mode, pe_data_mt->opt0.trim_qual);
            ct = clock();
        }
        {
        /*LPVOID *lp_sa[2];
        lp_sa[0] = (LPVOID*) ( (uint32_t*)pe_data_mt->lp_sa[0]+pe_data_mt->sai_i[0][b] );
        lp_sa[1] = (LPVOID*) ( (uint32_t*)pe_data_mt->lp_sa[1]+pe_data_mt->sai_i[1][b] );

        fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");
        cnt_chg = bwa_cal_pac_pos_pe_mt(pe_data_mt->bns, pe_data_mt->prefix, pe_data_mt->bwt, n_seqs, seqs,
            lp_sa, &ii, pe_data_mt->popt, &pe_data_mt->opt0, &last_ii);
        fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - ct) / CLOCKS_PER_SEC); ct = clock();
        fprintf(stderr, "[bwa_sai2sam_pe_core] changing coordinates of %d alignments.\n", cnt_chg);*/
        }//work1 ends

        //Signal_a next, rotate
        signal_t = (t+pe_data_mt->n_threads+1)%pe_data_mt->n_threads;
        fprintf(stderr, "T[%u] going to signal_a T[%u]\n", t, signal_t);

        SetEvent(pe_eid_a[signal_t]);

        //Outside serialized area, do concurrent works 1
        if (seqs[0] != 0)
        {
            LPVOID *lp_sa[2];
            lp_sa[0] = (LPVOID*) ( (uint32_t*)pe_data_mt->lp_sa[0]+pe_data_mt->sai_i[0][b] );
            lp_sa[1] = (LPVOID*) ( (uint32_t*)pe_data_mt->lp_sa[1]+pe_data_mt->sai_i[1][b] );

            fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");
            cnt_chg = bwa_cal_pac_pos_pe_mt(pe_data_mt->bns, pe_data_mt->prefix, pe_data_mt->bwt, n_seqs, seqs,
                lp_sa, &ii, pe_data_mt->popt, &pe_data_mt->opt0, &last_ii);
            fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - ct) / CLOCKS_PER_SEC); ct = clock();
            fprintf(stderr, "[bwa_sai2sam_pe_core] changing coordinates of %d alignments.\n", cnt_chg);
        
            fprintf(stderr, "[bwa_sai2sam_pe_core] align unmapped mate...\n");
            pacseq = bwa_paired_sw(pe_data_mt->bns, pe_data_mt->pac, n_seqs, seqs, pe_data_mt->popt, &ii);
            fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); ct = clock();

            fprintf(stderr, "[bwa_sai2sam_pe_core] refine gapped alignments... ");
            for (i = 0; i < 2; ++i)
                bwa_refine_gapped(pe_data_mt->bns, n_seqs, seqs[i], pacseq, pe_data_mt->ntbns);
            fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); ct = clock();
            if (pe_data_mt->pac == 0) free(pacseq);
        }
        //concurrent area 1 ends

        fprintf(stderr, "T[%u] going to wait on event_b\n", t);


        //wait for event b, except t0&b0, i.e. the first batch
        if (t != 0 || b != 0)
            WaitForSingleObject(pe_eid_b[t], INFINITE);

        //works2, serialized area
        /*fprintf(stderr, "T[%u] working_2 on batch [%u]\n", t, b);
        for(i=0;i<2;i++){
            fprintf(stderr, "[%u-%u]", t, b);
        }*/

        if (seqs[0] != 0)
        {
            fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
            for (i = 0; i < n_seqs; ++i) {
                bwa_seq_t *p[2];
                p[0] = seqs[0] + i; p[1] = seqs[1] + i;
                if (p[0]->bc[0] || p[1]->bc[0]) {
                    strcat(p[0]->bc, p[1]->bc);
                    strcpy(p[1]->bc, p[0]->bc);
                }
                bwa_print_sam1(pe_data_mt->bns, p[0], p[1], pe_data_mt->opt0.mode, pe_data_mt->opt0.max_top2);
                bwa_print_sam1(pe_data_mt->bns, p[1], p[0], pe_data_mt->opt0.mode, pe_data_mt->opt0.max_top2);
            }
            fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); ct = clock();

            tot_seqs += n_seqs; //moved here to stay close with print(), otherwise won't be accurate
            fprintf(stderr, "[bwa_sai2sam_pe_core] %d sequences have been processed.\n", tot_seqs); //moved here to stay close with add
        }
        //work2 ends

        //Signal_b next, rotate
        fprintf(stderr, "T[%u] going to signal_b T[%u]\n", t, signal_t);

        SetEvent(pe_eid_b[signal_t]);


        //Outside serialized area, do concurrent works 2
        for (i = 0; i < 2; ++i)
            bwa_free_read_seq(n_seqs, seqs[i]);
        

        b += pe_data_mt->n_threads; //next batch to process
        if (b > pe_data_mt->batch-1)
        {
            fprintf(stderr, "T[%u] exit\n", t);
            break; //no more batch to process, exit loop, thread exit
        }
        else
        {
            fprintf(stderr, "T[%u] going to loop over\n", t);
        }            
    } //end of while loop

    return 0;
}

typedef struct {
    kvec_t(bwt_aln1_t) aln;
} aln_buf_t;

int bwa_cal_pac_pos_pe_mt(const bntseq_t *bns, const char *prefix, bwt_t *const _bwt, int n_seqs, bwa_seq_t *seqs[2], LPVOID *lp_sa[2], isize_info_t *ii,
                       const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii)
{
    int i, j, cnt_chg = 0;
    char str[1024];
    bwt_t *bwt;
    pe_data_t *d;
    aln_buf_t *buf[2];
    int sa_offset[2] = {0,0}; //to keep SAIs offset

    int const tsize = sizeof(bwt_aln1_t);

    d = (pe_data_t*)calloc(1, sizeof(pe_data_t));
    buf[0] = (aln_buf_t*)calloc(n_seqs, sizeof(aln_buf_t));
    buf[1] = (aln_buf_t*)calloc(n_seqs, sizeof(aln_buf_t));

    if (_bwt == 0) { // load forward SA
        strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
        strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
    } else bwt = _bwt;

    // SE
    for (i = 0; i != n_seqs; ++i) {
        bwa_seq_t *p[2];
        for (j = 0; j < 2; ++j) {
            int n_aln;
            p[j] = seqs[j] + i;
            p[j]->n_multi = 0;
            p[j]->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);
            n_aln = *((uint32_t*)lp_sa[j]+sa_offset[j]++); //new, take n_aln, then move 1 step

            if (n_aln > kv_max(d->aln[j]))
                kv_resize(bwt_aln1_t, d->aln[j], n_aln);
            d->aln[j].n = n_aln; 
            
            memcpy(d->aln[j].a, (uint32_t*)lp_sa[j]+sa_offset[j], n_aln * tsize); //new, take alns
            sa_offset[j] += n_aln * tsize/4; // sizeof(uint32_t) is 4

            kv_copy(bwt_aln1_t, buf[j][i].aln, d->aln[j]); // backup d->aln[j]
            // generate SE alignment and mapping quality
            bwa_aln2seq(n_aln, d->aln[j].a, p[j]);
            if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {
                int strand;
                int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
                p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
                p[j]->pos = bwa_sa2pos(bns, bwt, p[j]->sa, p[j]->len, &strand);
                p[j]->strand = strand;
            }
        }
    }

    // infer isize
    infer_isize(n_seqs, seqs, ii, opt->ap_prior, bwt->seq_len/2);
    if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
    if (opt->force_isize) {
        fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __FUNCTION__);
        ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
    }

    // PE
    for (i = 0; i != n_seqs; ++i) {
        bwa_seq_t *p[2];
        for (j = 0; j < 2; ++j) {
            p[j] = seqs[j] + i;
            kv_copy(bwt_aln1_t, d->aln[j], buf[j][i].aln);
        }
        if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
            && (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
        { // only when both ends mapped
            khint_t iter;
            b128_t x;
            int j, k;
            long long n_occ[2];
            for (j = 0; j < 2; ++j) {
                n_occ[j] = 0;
                for (k = 0; k < d->aln[j].n; ++k)
                    n_occ[j] += d->aln[j].a[k].l - d->aln[j].a[k].k + 1;
            }
            if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;
            d->arr.n = 0;
            for (j = 0; j < 2; ++j) {
                for (k = 0; k < d->aln[j].n; ++k) {
                    bwt_aln1_t *r = d->aln[j].a + k;
                    bwtint_t l;
                    if (0 && r->l - r->k + 1 >= MIN_HASH_WIDTH) { // then check hash table
                        b128_t key;
                        int ret;
                        key.x = r->k; key.y = r->l;
                        iter = kh_put(b128, g_hash, key, &ret);
                        if (ret) { // not in the hash table; ret must equal 1 as we never remove elements
                            poslist_t *z = &kh_val(g_hash, iter);
                            z->n = r->l - r->k + 1;
                            z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
                            for (l = r->k; l <= r->l; ++l) {
                                int strand;
                                z->a[l - r->k] = bwa_sa2pos(bns, bwt, l, p[j]->len, &strand)<<1;
                                z->a[l - r->k] |= strand;
                            }
                        }
                        for (l = 0; l < kh_val(g_hash, iter).n; ++l) {
                            x.x = kh_val(g_hash, iter).a[l]>>1;
                            x.y = k<<2 | (kh_val(g_hash, iter).a[l]&1)<<1 | j;
                            kv_push(b128_t, d->arr, x);
                        }
                    } else { // then calculate on the fly
                        for (l = r->k; l <= r->l; ++l) {
                            int strand;
                            x.x = bwa_sa2pos(bns, bwt, l, p[j]->len, &strand);
                            x.y = k<<2 | strand<<1 | j;
                            kv_push(b128_t, d->arr, x);
                        }
                    }
                }
            }
            cnt_chg += pairing(p, d, opt, gopt->s_mm, ii);
        }

        if (opt->N_multi || opt->n_multi) {
            for (j = 0; j < 2; ++j) {
                if (p[j]->type != BWA_TYPE_NO_MATCH) {
                    int k, n_multi;
                    if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
                        bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0, p[j]->c1+p[j]->c2-1 > opt->N_multi? opt->n_multi : opt->N_multi);
                    } else bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0, opt->n_multi);
                    for (k = 0, n_multi = 0; k < p[j]->n_multi; ++k) {
                        int strand;
                        bwt_multi1_t *q = p[j]->multi + k;
                        q->pos = bwa_sa2pos(bns, bwt, q->pos, p[j]->len, &strand);
                        q->strand = strand;
                        if (q->pos != p[j]->pos)
                            p[j]->multi[n_multi++] = *q;
                    }
                    p[j]->n_multi = n_multi;
                }
            }
        }
    }

    // free
    for (i = 0; i < n_seqs; ++i) {
        kv_destroy(buf[0][i].aln);
        kv_destroy(buf[1][i].aln);
    }
    free(buf[0]); free(buf[1]);
    if (_bwt == 0) bwt_destroy(bwt);
    kv_destroy(d->arr);
    kv_destroy(d->pos[0]); kv_destroy(d->pos[1]);
    kv_destroy(d->aln[0]); kv_destroy(d->aln[1]);
    free(d);
    return cnt_chg;
}

//Modified from bwtio.c
bwt_t *bwt_restore_bwt_mmf(char *fn)
{
    bwt_t *bwt;
    LARGE_INTEGER  lpFileSize;
    LPVOID lpView; BOOL b;
    int ts = sizeof(bwtint_t);
    uint64_t i = 0, j = 0;

    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    lpView = mmf_open_ex(fn, &lpFileSize);
    bwt->bwt_size = (lpFileSize.QuadPart - ts * 5) >> 2;
    bwt->bwt = (uint32_t*)lpView + ts * 5 / 4; //uint32_t is 4 bytes
    memcpy(&bwt->primary, (uint32_t*)lpView, ts);
    memcpy(bwt->L2+1, (uint32_t*)lpView+2, ts * 4); //bwtint_t is 8 bytes, step is 8/4*1
    bwt->seq_len = bwt->L2[4];
    bwt_gen_cnt_table(bwt);

    for(i = 0; i < lpFileSize.QuadPart/4; i++)
        j += bwt->bwt[i] & 0x3; //touch every byte to commit into memory
    i = j; //useless, just trying to convince compiler don't ignore the loop

    b = VirtualLock(lpView, lpFileSize.QuadPart);
    if (!b) fprintf(stderr, "Failed to lock for bwt\n");

    return bwt;
}

//Modified from bwtio.c
void bwt_restore_sa_mmf(char *fn, bwt_t *bwt)
{
    bwtint_t primary;
    LARGE_INTEGER  lpFileSize;
    LPVOID lpView; BOOL b;
    int ts = sizeof(bwtint_t);
    uint64_t i = 0, j = 0;

    lpView = mmf_open_ex(fn, &lpFileSize);

    memcpy(&primary, (uint64_t*)lpView, ts);
    xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
    memcpy(&bwt->sa_intv, (uint64_t*)lpView+5, ts); //skip 4, plus the first 1
    memcpy(&primary, (uint64_t*)lpView+6, ts); // 1 after above, primary as variable is just re-used here
                                               // as a placeholder, it actually is holding seq_len
    xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

    bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
    bwt->sa = (bwtint_t*)lpView + 6; //SA starts from 7, but that's for sa[1],
                                     // my workaround is use 6, that basically means sa[0] = seq_len

    //bwt->sa[0] = -1; //my mmf is readonly, I can't do this setting
    for(i = 0; i < bwt->n_sa; i++)
        j += bwt->sa[i] & 0x3; //touch every byte to commit into memory
    i = j;

    b = VirtualLock(lpView, lpFileSize.QuadPart);
    if (!b) fprintf(stderr, "Failed to lock for sa\n");
}

#endif