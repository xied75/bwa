/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "utils.h"
#include "bwt.h"

static const uint64_t occ_mask[32] = {
	0xc0000000ul, 	0xf0000000ul, 	0xfc000000ul,
	0xff000000ul, 	0xffc00000ul, 	0xfff00000ul,
	0xfffc0000ul, 	0xffff0000ul, 	0xffffc000ul,
	0xfffff000ul, 	0xfffffc00ul, 	0xffffff00ul,
	0xffffffc0ul, 	0xfffffff0ul, 	0xfffffffcul,
	0xfffffffful, 	0xc0000000fffffffful, 	0xf0000000fffffffful,
	0xfc000000fffffffful, 	0xff000000fffffffful, 	0xffc00000fffffffful,
	0xfff00000fffffffful, 	0xfffc0000fffffffful, 	0xffff0000fffffffful,
	0xffffc000fffffffful, 	0xfffff000fffffffful, 	0xfffffc00fffffffful,
	0xffffff00fffffffful, 	0xffffffc0fffffffful, 	0xfffffff0fffffffful,
	0xfffffffcfffffffful, 0xfffffffffffffffful
};
static const uint64_t n_mask[5] = { 0xfffffffffffffffful, 0xaaaaaaaaaaaaaaaaul, 
		0x5555555555555555ul, 0x0ul, 0xfffffffffffffffful };

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

static inline int __occ_aux(uint64_t y)
{
	// reduce nucleotide counting to bits counting
	y = (y >> 1) & y;
	// count the number of 1s in y
	y = (y & 0x1111111111111111ul) + (y >> 2 & 0x1111111111111111ul);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0ful) * 0x101010101010101ul >> 56;
	// this, instead of return statement, also works:
	//y = (y + (y >> 4u)) & 0x707070707070707ul; // at most 4 per f
	//y = (y + (y >> 8u)) & 0xf000f000f000ful; // at most 8 per f
	//y = (y + (y >> 16u)) & 0x1f0000001ful; // at most 16 per 1f
	//return y + (y >> 32u);
}

#define __occ_aux_p(z) ({ 						\
	z = (z >> 1) & z;						\
	(z & 0x1111111111111111ul) + (z >> 2 & 0x1111111111111111ul);	\
})

#define __occ_aux_p2(y) ({						\
	y = y + (y >> 4u);						\
	(y & 0xf0f0f0f0f0f0f0ful) * 0x101010101010101ul >> 56;		\
/*	y = (y + (y >> 4u)) & 0x0f0f0f0f0f0f0f0ful;			\
	y = (y + (y >> 8u)) & 0x1f001f001f001ful;			\
	y = (y + (y >> 16u)) & 0x3f0000003ful;				\
	y + (y >> 32u);					*/		\
})
//	y = y + (y >> 4u);						\
//	(y & 0xf0f0f0f0f0f0f0ful) * 0x101010101010101ul >> 56;		\

#define bwt_occ_pn(z, y, l, n, trdp)				\
	switch (l) {						\
		case 3: (z) = trdp;				\
			(z) = __occ_aux_p(z);			\
			n += __occ_aux_p2(z);			\
		case 2: (z) = trdp;				\
			(y) += __occ_aux_p(z);			\
		case 1: (z) = trdp;				\
			(y) += __occ_aux_p(z);			\
	}							\

#define bwt_occ_p(z, y, l, trdp)				\
	switch (l) {						\
		case 3: (z) = trdp;				\
			(y) += __occ_aux_p(z);			\
		case 2: (z) = trdp;				\
			(y) += __occ_aux_p(z);			\
		case 1: (z) = trdp;				\
			(y) += __occ_aux_p(z);			\
	}

static inline bwtint_t bwt_occ_0(bwtint_t l, const uint64_t *p)
{
	uint64_t z, y = 0ul;
	switch (l) {
		case 3: z = -*(--p) - 1ul;
			y += __occ_aux_p(z);
		case 2: z = -*(--p) - 1ul;
			y += __occ_aux_p(z);
		default: z = -*(--p) - 1ul;
			y += __occ_aux_p(z);
	}
	return __occ_aux_p2(y);
}

static inline bwtint_t bwt_occ_x(const uint64_t x, bwtint_t l, const uint64_t *p)
{
	uint64_t z, y = 0ul;
	switch (l) {
		case 3: z = *(--p) ^ x;
			y += __occ_aux_p(z);
		case 2: z = *(--p) ^ x;
			y += __occ_aux_p(z);
		default: z = *(--p) ^ x;
			y += __occ_aux_p(z);
	}
	return __occ_aux_p2(y);
}

static inline bwtint_t bwt_occ_3(bwtint_t l, const uint64_t *p)
{
	uint64_t z, y = 0ul;
	switch (l) {
		case 3: z = *(--p);
			y += __occ_aux_p(z);
		case 2: z = *(--p);
			y += __occ_aux_p(z);
		default: z = *(--p);
			y += __occ_aux_p(z);
	}
	return __occ_aux_p2(y);
}

static inline bwtint_t bwt_occ(const uint64_t *p, const bwtint_t ko, bwtint_t k, const ubyte_t c)
{
	uint64_t z, y = 0ul;
	bwtint_t n, l;
	p += ko * 6;
	n = ((uint32_t *)p)[c];
	k -= ko * OCC_INTERVAL;
	l = k >> 5;
	p += 2; // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	switch (c) {
		case 0: k &= 31;
			bwt_occ_pn(z, y, l, n, -*(p++) - 1ul)
			z = -(*p & occ_mask[k]) - 1ul;
			n -= (k^31);
			break;
		case 3: bwt_occ_pn(z, y, l, n, *(p++))
			z = *p & occ_mask[k&31];
			break;
		default: bwt_occ_pn(z, y, l, n, *(p++) ^ n_mask[c])
			z = (*p & occ_mask[k&31]) ^ n_mask[c];
	}
	y += __occ_aux_p(z);
	n += __occ_aux_p2(y);
	return n;
}

static inline bwtint_t cal_isa(const bwt_t *bwt, bwtint_t isa)
{
	bwtint_t c, _isa, isa_o;
	if (isa != bwt->primary) {
		_isa = (isa < bwt->primary) ? isa : isa - 1;
		isa_o = _isa/OCC_INTERVAL;
		c = bwt_B0(bwt, _isa, isa_o);
		if (isa < bwt->seq_len) {
			isa = bwt->L2[c] + bwt_occ((uint64_t *)bwt->bwt, isa_o, _isa, c);
		} else {
			isa = (isa == bwt->seq_len ? bwt->L2[c+1] : bwt->L2[c]);
		}
	} else {
		isa = 0;
	}
	return isa;
}

// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa

	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	// calculate SA value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0) bwt->sa[isa/intv] = sa;
		--sa;
		isa = cal_isa(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t m, add, sa = 0;
	m = bwt->sa_intv - 1;
	if ((m+1) & m) { // not power of 2 before decrement
		add = m;
		m |= m>>1;
		m |= m>>2;
		m |= m>>4;
		m |= m>>8;
		m |= m>>16;
		add ^= m;
	} else {
		add = 0;
	}
	while (((k + add) & m)) {
		++sa;
		k = cal_isa(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
inline bwtint_t bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t *l, ubyte_t c)
{
	bwtint_t i, m, n, kr;
	const uint64_t *p;
	n = k - 1;
	m = (n >= bwt->primary)? n-1 : n;
	i = (*l >= bwt->primary)? *l-1 : *l;
	n = m/OCC_INTERVAL;
	*l = i/OCC_INTERVAL;
	kr = bwt->L2[c];
	if (*l == n && k != 0) {
		uint64_t z, y;
		p = (const uint64_t *)bwt->bwt + n * 6;
		k = m & 31;
		m = ((m - (n * OCC_INTERVAL)) >> 5); //k-len
		n = i & 31;
		i = ((i - (*l * OCC_INTERVAL)) >> 5); //total len
		kr += ((uint32_t *)p)[c];
		p += 2 + i; // jump to the end of the last BWT cell
		i = i - m; //total len - k len = l-len
		if (c != 0) {
			z = (*p & occ_mask[n]) ^ n_mask[c];
			y = __occ_aux_p(z);
			bwt_occ_pn(z, y, i, *l, *(--p) ^ n_mask[c])
			*l = __occ_aux_p2(y);

			z = (*p & occ_mask[k]) ^ n_mask[c];
			y = __occ_aux_p(z);
			k = __occ_aux_p2(y) + 1;
			if (k > *l) {
				*l = bwt->seq_len;
				return 0;
			}
			y = 0ul;
			bwt_occ_p(z, y, m, *(--p) ^ n_mask[c])
			kr += __occ_aux_p2(y);
		} else {
			z = -(*p & occ_mask[n]) -1ul;
			y = __occ_aux_p(z);
			bwt_occ_pn(z, y, i, *l, -*(--p) -1ul)
			*l = __occ_aux_p2(y) - (n^31);

			z = -(*p & occ_mask[k]) -1ul;
			y = __occ_aux_p(z);
			k = __occ_aux_p2(y) - (k^31) + 1;
			if (k > *l) {
				*l = bwt->seq_len;
				return 0;
			}
			y = 0ul;
			bwt_occ_p(z, y, m, -*(--p) -1ul)
			kr += __occ_aux_p2(y);
		}
		*l += kr;
		kr += k;
	} else {
		if (k != 0) {
			*l = kr + bwt_occ((const uint64_t *)bwt->bwt, *l, i, c);
			kr += bwt_occ((const uint64_t *)bwt->bwt, n, m, c);
			if (++kr > *l) {
				*l = bwt->seq_len;
				kr = 0;
			}
		} else {
			*l = bwt->L2[c+1];
			++kr;
		}
	}
	return kr;
}

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t l, j, x;
	uint32_t *p;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	if (k >= bwt->primary) --k; // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 16);
	p += 4;
	j = k >> 4 << 4;
	for (l = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; l < j; l += 16, ++p)
		x += __occ_aux4(bwt, *p);
	x += __occ_aux4(bwt, *p & occ_mask[k&15]) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	if (k == l) {
		bwt_occ4(bwt, k, cntk);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		return;
	}
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t i, j, x, y;
		uint32_t *p;
		int cl[4];
		if (k >= bwt->primary) --k; // because $ is not in bwt
		if (l >= bwt->primary) --l;
		cl[0] = cl[1] = cl[2] = cl[3] = 0;
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += 4;
		// prepare cntk[]
		j = k >> 4 << 4;
		for (i = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; i < j; i += 16, ++p)
			x += __occ_aux4(bwt, *p);
		y = x;
		x += __occ_aux4(bwt, *p & occ_mask[k&15]) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		j = l >> 4 << 4;
		for (; i < j; i += 16, ++p)
			y += __occ_aux4(bwt, *p);
		y += __occ_aux4(bwt, *p & occ_mask[l&15]) - (~l&15);
		memcpy(cntl, cntk, 16);
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3 || !(k = bwt_2occ(bwt, k, &l, c)))
			return 0; // no match
	}
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3 || !(k = bwt_2occ(bwt, k, &l, c)))
			return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}
