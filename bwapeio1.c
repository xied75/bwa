#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

extern int num_sampe_threads;
extern int async_read_seq;

#ifdef HAVE_PTHREAD
pthread_mutex_t io_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_t read_tid = 0;
#endif // HAVE_PTHREAD

static bwa_seqio_t *ks1_copy = NULL;
static bwa_seqio_t *ks2_copy = NULL;
static int *n_read_addr = 0;
static int mode1_copy = 0;
static int mode2_copy = 0;
static int trim1_copy = 0;
static int trim2_copy = 0;
static bwa_seq_t **seq1_addr = NULL;
static bwa_seq_t **seq2_addr = NULL;


// -------------------

static void thr_bwa_read_seq2_tpx(long n_needed)
{
	int n1, n2;

	*seq1_addr = bwa_read_seq(ks1_copy, (int)n_needed, &n1, mode1_copy, trim1_copy);

	if(*seq1_addr == NULL){
		return;
	}

	*seq2_addr = bwa_read_seq(ks2_copy, (int)n_needed, &n2, mode2_copy, trim2_copy);

	if(n1 != n2){
		fflush(NULL);
		fprintf(stderr,"[tpx_bwa_read_seq2] n1 (%d) != n2 (%d)\n",n1,n2);
		exit(1);
	}

	*n_read_addr = n1;

	return;
}

// -------------------

void bwa_read_seq2_wait_tpx(void)
{
#ifdef HAVE_PTHREAD
        if( (async_read_seq) && (read_tid != 0) ){
                pthread_join(read_tid, NULL);
                read_tid = 0;
        }
#endif // HAVE_PTHREAD

        return;
}

// -------------------

void bwa_read_seq2_tpx(bwa_seqio_t *ks1, bwa_seqio_t *ks2, int n_needed, int *n,
                       int mode1, int mode2, int trim_qual1, int trim_qual2,
                       bwa_seq_t **seq1, bwa_seq_t **seq2)
{
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

#ifdef _DEBUG
	pthread_mutex_lock(&io_lock);
	fprintf(stderr,"bwapeio1: n_needed = %d\n",n_needed);
	pthread_mutex_unlock(&io_lock);
#endif

#ifdef _TIMING
	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
#endif

#ifdef HAVE_PTHREAD
	if( (async_read_seq) && (num_sampe_threads > 1) ){

		if(read_tid != 0){
			pthread_join(read_tid, NULL);
			read_tid = 0;
		}

		ks1_copy = ks1;
		ks2_copy = ks2;
		n_read_addr = n;
		mode1_copy = mode1;
		mode2_copy = mode2;
		trim1_copy = trim_qual1;
		trim2_copy = trim_qual2;
		seq1_addr = seq1;
		seq2_addr = seq2;

		int srtn;
		srtn = pthread_create(&read_tid,NULL,(void *(*)(void *))thr_bwa_read_seq2_tpx,(void *)(long)n_needed);
		if(srtn != 0){
			fprintf(stderr,"[bwa_read_seq2_tpx] pthread_create thr_bwa_read_seq2_tpx error %d\n",srtn);
			exit(1);
		}

	}else{

		ks1_copy = ks1;
		ks2_copy = ks2;
		n_read_addr = n;
		mode1_copy = mode1;
		mode2_copy = mode2;
		trim1_copy = trim_qual1;
		trim2_copy = trim_qual2;
		seq1_addr = seq1;
		seq2_addr = seq2;

		thr_bwa_read_seq2_tpx(n_needed);

	}
#else // HAVE_PTHREAD
	ks1_copy = ks1;
	ks2_copy = ks2;
	n_read_addr = n;
	mode1_copy = mode1;
	mode2_copy = mode2;
	trim1_copy = trim_qual1;
	trim2_copy = trim_qual2;
	seq1_addr = seq1;
	seq2_addr = seq2;

	thr_bwa_read_seq2_tpx(n_needed);
#endif // HAVE_PTHREAD

#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
	pthread_mutex_lock(&io_lock);
# endif // HAVE_PTHREAD
	fprintf(stderr,"bwapeio1 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&io_lock);
# endif // HAVE_PTHREAD
#endif

	return;
}
