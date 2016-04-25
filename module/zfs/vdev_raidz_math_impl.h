/*
 * CDDL HEADER START
 *
 * The contents of this file are subject to the terms of the
 * Common Development and Distribution License (the "License").
 * You may not use this file except in compliance with the License.
 *
 * You can obtain a copy of the license at usr/src/OPENSOLARIS.LICENSE
 * or http://www.opensolaris.org/os/licensing.
 * See the License for the specific language governing permissions
 * and limitations under the License.
 *
 * When distributing Covered Code, include this CDDL HEADER in each
 * file and include the License file at usr/src/OPENSOLARIS.LICENSE.
 * If applicable, add the following below this CDDL HEADER, with the
 * fields enclosed by brackets "[]" replaced with your own identifying
 * information: Portions Copyright [yyyy] [name of copyright owner]
 *
 * CDDL HEADER END
 */
/*
 * Copyright (C) 2016 Gvozden Nešković. All rights reserved.
 */

#ifndef _VDEV_RAIDZ_MATH_IMPL_H
#define	_VDEV_RAIDZ_MATH_IMPL_H

#include <sys/types.h>

#define	RAIDZ_BUG()	ASSERT(0)

#define	raidz_inline inline __attribute__((always_inline))

#define	COL_OFF(col, off)	((v_t *)(((char *)(col)->rc_data) + (off)))

/*
 * PARITY CALCULATION
 * An optimized function is called for a full length of data columns
 * If RAIDZ map contains remainder columns (shorter columns) the same function
 * is called for reminder of full columns.
 *
 * GEN_[P|PQ|PQR]_BLOCK() functions are designed to be efficiently in-lined by
 * the compiler. This removes a lot of conditionals from the inside loop which
 * makes the code faster, especially for vectorized code.
 * They are also highly parametrized, allowing for each implementation to define
 * most optimal stride, and register allocation.
 */

static raidz_inline void
GEN_P_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const boolean_t big)
{
	int c;
	size_t ioff;
	const size_t ncols = big ? raidz_nbigcols(rm) : raidz_ncols(rm);
	raidz_col_t * const pcol = &(rm->rm_col[CODE_P]);
	raidz_col_t *col;

	GEN_P_DEFINE();

	for (ioff = off; ioff < end; ioff += (GEN_P_STRIDE * sizeof (v_t))) {
		LOAD(COL_OFF(&(rm->rm_col[1]), ioff), GEN_P_P);
		for (c = 2; c < ncols; c++) {
			col = &rm->rm_col[c];
			XOR_ACC(COL_OFF(col, ioff), GEN_P_P);
		}
		STORE(COL_OFF(pcol, ioff), GEN_P_P);
	}
}

/*
 * Generate P parity (RAIDZ1)
 *
 * @rm	RAIDZ map
 */
static raidz_inline void
raidz_generate_p_impl(raidz_map_t * const rm)
{
	const size_t ncols = raidz_ncols(rm);
	const size_t psize = rm->rm_col[CODE_P].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;

	raidz_math_begin();

	/* lastcsize */
	GEN_P_BLOCK(rm, 0, lastcsize, B_FALSE);

	/* fullcols */
	if (lastcsize < psize)
		GEN_P_BLOCK(rm, lastcsize, psize, B_TRUE);

	raidz_math_end();
}

static raidz_inline void
GEN_PQ_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const boolean_t big)
{
	int c;
	size_t ioff;
	const size_t firstdc = raidz_parity(rm);
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = big ? raidz_nbigcols(rm) : ncols;
	raidz_col_t * const pcol = &rm->rm_col[CODE_P];
	raidz_col_t * const qcol = &rm->rm_col[CODE_Q];
	raidz_col_t *col;

	GEN_PQ_DEFINE();

	MUL2_SETUP();

	for (ioff = off; ioff < end; ioff += (GEN_PQ_STRIDE * sizeof (v_t))) {
		LOAD(COL_OFF(&rm->rm_col[firstdc], ioff), GEN_PQ_P);
		COPY(GEN_PQ_P, GEN_PQ_Q);
		for (c = firstdc+1; c < nbigcols; c++) {
			col = &rm->rm_col[c];
			LOAD(COL_OFF(col, ioff), GEN_PQ_D);
			MUL2(GEN_PQ_Q);
			XOR(GEN_PQ_D, GEN_PQ_P);
			XOR(GEN_PQ_D, GEN_PQ_Q);
		}
		STORE(COL_OFF(pcol, ioff), GEN_PQ_P);
		if (big) {
			for (; c < ncols; c++)
				MUL2(GEN_PQ_Q);
		}
		STORE(COL_OFF(qcol, ioff), GEN_PQ_Q);
	}
}

/*
 * Generate PQ parity (RAIDZ2)
 *
 * @rm	RAIDZ map
 */
static raidz_inline void
raidz_generate_pq_impl(raidz_map_t * const rm)
{
	const size_t ncols = raidz_ncols(rm);
	const size_t psize = rm->rm_col[CODE_P].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;

	raidz_math_begin();

	/* lastcsize */
	GEN_PQ_BLOCK(rm, 0, lastcsize, B_FALSE);

	/* fullcols */
	if (lastcsize < psize) {
		GEN_PQ_BLOCK(rm, lastcsize, psize, B_TRUE);
	}

	raidz_math_end();
}


static raidz_inline void
GEN_PQR_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const boolean_t big)
{
	int c;
	size_t ioff;
	const size_t firstdc = raidz_parity(rm);
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = big ? raidz_nbigcols(rm) : ncols;
	raidz_col_t *col;
	raidz_col_t * const pcol = &rm->rm_col[CODE_P];
	raidz_col_t * const qcol = &rm->rm_col[CODE_Q];
	raidz_col_t * const rcol = &rm->rm_col[CODE_R];

	GEN_PQR_DEFINE();

	MUL2_SETUP();

	for (ioff = off; ioff < end; ioff += (GEN_PQR_STRIDE * sizeof (v_t))) {
		LOAD(COL_OFF(&rm->rm_col[firstdc], ioff), GEN_PQR_P);
		COPY(GEN_PQR_P, GEN_PQR_Q);
		COPY(GEN_PQR_P, GEN_PQR_R);
		for (c = firstdc+1; c < nbigcols; c++) {
			col = &rm->rm_col[c];
			LOAD(COL_OFF(col, ioff), GEN_PQR_D);
			MUL2(GEN_PQR_Q);
			MUL4(GEN_PQR_R);
			XOR(GEN_PQR_D, GEN_PQR_P);
			XOR(GEN_PQR_D, GEN_PQR_Q);
			XOR(GEN_PQR_D, GEN_PQR_R);
		}
		STORE(COL_OFF(pcol, ioff), GEN_PQR_P);

		if (big) {
			for (; c < ncols; c++) {
				MUL2(GEN_PQR_Q);
				MUL4(GEN_PQR_R);
			}
		}

		STORE(COL_OFF(qcol, ioff), GEN_PQR_Q);
		STORE(COL_OFF(rcol, ioff), GEN_PQR_R);
	}
}


/*
 * Generate PQR parity (RAIDZ3)
 *
 * @rm	RAIDZ map
 */
static raidz_inline void
raidz_generate_pqr_impl(raidz_map_t * const rm)
{
	const size_t ncols = raidz_ncols(rm);
	const size_t psize = rm->rm_col[CODE_P].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;

	raidz_math_begin();

	/* lastcsize */
	GEN_PQR_BLOCK(rm, 0, lastcsize, B_FALSE);

	/* fullcols */
	if (lastcsize < psize)
		GEN_PQR_BLOCK(rm, lastcsize, psize, B_TRUE);

	raidz_math_end();
}


/*
 * DATA RECONSTRUCTION
 *
 * Data reconstruction process consists of two phases:
 * 	- Syndrome calculation
 * 	- Data reconstruction
 *
 * Syndrome is calculated by generating parity using available data columns
 * and zeros in places of erasure. Finally, Existing parity is added to
 * corresponding syndrome value to obtain the [P|Q|R]syn values from equation:
 * 	P = Psyn + Dx + Dy + Dz
 * 	Q = Qsyn + 2^x * Dx + 2^y * Dy + 2^z * Dz
 * 	R = Rsyn + 4^x * Dx + 4^y * Dy + 4^z * Dz
 *
 * For the data reconstruction phase, the corresponding equations are solved
 * for Dx, Dy, Dz (missing data). This generally involves multiplying known
 * symbols by an coefficient and adding them together. The multiplication
 * constant coefficients are calculated ahead of the operation in
 * raidz_init_rec_coeff().
 *
 * IMPLEMENTATION NOTE:
 * For linear zio buffers, where all parts of the buffer are accessible,
 * syndrome calculation and data reconstruction phases can be performed in
 * the single pass over the zio buffer.
 */


/*
 * Function calculates multiplication constants for specified reconstruction.
 * Coefficients depend on RAIDZ geometry, indexes of failed child vdevs, and
 * used parity columns for reconstruction.
 *
 * @rm			RAIDZ map
 * @tgtidx		array of missing data indexes
 * @raidz_rec_op	reconstruct operation (which parity is used)
 * @coeff		output array of coefficients. Array must be user
 *         		provided and must hold minimum MUL_CNT values
 */
static raidz_inline void
raidz_init_rec_coeff(const raidz_map_t * const rm, const int *tgtidx,
    const enum raidz_rec_op op, unsigned *coeff)
{
	const unsigned ncols = raidz_ncols(rm);
	const unsigned x = tgtidx[TARGET_X];
	const unsigned y = tgtidx[TARGET_Y];
	const unsigned z = tgtidx[TARGET_Z];
	unsigned a, b, e, denom, x_d, y_d;

	switch (op) {
		case RAIDZ_REC_P:
			break;
		case RAIDZ_REC_Q:
			coeff[MUL_Q_X] = fix_mul_exp(255 - (ncols - x - 1));
			break;
		case RAIDZ_REC_R:
			coeff[MUL_R_X] = fix_mul_exp(255 - 2 * (ncols - x - 1));
			break;
		case RAIDZ_REC_PQ:
			a = vdev_raidz_pow2[255 + x - y];
			b = vdev_raidz_pow2[255 - (ncols - 1 - x)];
			e = 255 - vdev_raidz_log2[a ^ 0x01];

			coeff[MUL_PQ_X] = fix_mul_exp(
			    vdev_raidz_log2[vdev_raidz_exp2(a, e)]);
			coeff[MUL_PQ_Y] = fix_mul_exp(
			    vdev_raidz_log2[vdev_raidz_exp2(b, e)]);
			break;

		case RAIDZ_REC_PR:
			a = vdev_raidz_pow2[255 + 2 * x - 2 * y];
			b = vdev_raidz_pow2[255 - 2 * (ncols - 1 - x)];
			e = 255 - vdev_raidz_log2[a ^ 0x01];

			coeff[MUL_PR_X] = fix_mul_exp(
			    vdev_raidz_log2[vdev_raidz_exp2(a, e)]);
			coeff[MUL_PR_Y] = fix_mul_exp(
			    vdev_raidz_log2[vdev_raidz_exp2(b, e)]);
			break;

		case RAIDZ_REC_QR:
			denom = 255 - vdev_raidz_log2[
			    vdev_raidz_pow2[3 * ncols - 3 - x - 2 * y] ^
			    vdev_raidz_pow2[3 * ncols - 3 - 2 * x - y]
			];

			coeff[MUL_QR_XQ] = fix_mul_exp(ncols - 1 - y);
			coeff[MUL_QR_X]	= fix_mul_exp(ncols - 1 - y + denom);
			coeff[MUL_QR_YQ] = fix_mul_exp(ncols - 1 - x);
			coeff[MUL_QR_Y]	= fix_mul_exp(ncols - 1 - x + denom);
			break;

		case RAIDZ_REC_PQR:
			x_d = 255 - vdev_raidz_log2[
			    vdev_raidz_pow2[3 * ncols - 3 - 2 * x - y] ^
			    vdev_raidz_pow2[3 * ncols - 3 - x - 2 * y] ^
			    vdev_raidz_pow2[3 * ncols - 3 - 2 * x - z] ^
			    vdev_raidz_pow2[3 * ncols - 3 - x - 2 * z] ^
			    vdev_raidz_pow2[3 * ncols - 3 - 2 * y - z] ^
			    vdev_raidz_pow2[3 * ncols - 3 - y - 2 * z]
			];
			y_d = 255 - vdev_raidz_log2[
			    vdev_raidz_pow2[ncols - 1 - y] ^
			    vdev_raidz_pow2[ncols - 1 - z]
			];

			coeff[MUL_PQR_XP] = fix_mul_exp(vdev_raidz_log2[
			    vdev_raidz_pow2[3 * ncols - 3 - 2 * y - z] ^
			    vdev_raidz_pow2[3 * ncols - 3 - y - 2 * z]
			] + x_d);
			coeff[MUL_PQR_XQ] = fix_mul_exp(vdev_raidz_log2[
			    vdev_raidz_pow2[2 * ncols - 2 - 2 * y] ^
			    vdev_raidz_pow2[2 * ncols - 2 - 2 * z]
			] + x_d);
			coeff[MUL_PQR_XR] = fix_mul_exp(vdev_raidz_log2[
			    vdev_raidz_pow2[ncols - 1 - y] ^
			    vdev_raidz_pow2[ncols - 1 - z]
			] + x_d);
			coeff[MUL_PQR_YU] = fix_mul_exp(ncols - 1 - x);
			coeff[MUL_PQR_YP] = fix_mul_exp(ncols - 1 - z + y_d);
			coeff[MUL_PQR_YQ] = fix_mul_exp(y_d);
			break;
		default:
			RAIDZ_BUG();
			break;
	}
}

static raidz_inline void
REC_P_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const size_t ncols, const size_t nfullcols)
{
	int c;
	size_t ioff;
	const int x = tgtidx[TARGET_X];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const pcol = &rm->rm_col[CODE_P];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t *col;

	REC_P_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_P_STRIDE * sizeof (v_t))) {
		LOAD(COL_OFF(pcol, ioff), REC_P_X);
		for (c = firstdc; c < nfullcols; c++) {
			if (c != x) {
				col = &rm->rm_col[c];
				XOR_ACC(COL_OFF(col, ioff), REC_P_X);
			}
		}
		STORE(COL_OFF(xcol, ioff), REC_P_X);
	}
}

/*
 * Reconstruct single data column using P parity
 * @rec_method	REC_P_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */

static raidz_inline int
raidz_reconstruct_p_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t fullcsize = rm->rm_col[ncols-1].rc_size;

	raidz_math_begin();

	/* fullcsize */
	REC_P_BLOCK(rm, 0, fullcsize, tgtidx, ncols, ncols);

	/* fullcols */
	if (fullcsize < xsize)
		REC_P_BLOCK(rm, fullcsize, xsize, tgtidx, ncols, nbigcols);

	raidz_math_end();

	return (1 << CODE_P);
}


static raidz_inline void
REC_Q_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const unsigned *coeff, const size_t ncols,
    const size_t nbigcols)
{
	int c;
	size_t ioff = 0;
	const int x = tgtidx[TARGET_X];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const qcol = &rm->rm_col[CODE_Q];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t *col;

	REC_Q_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_Q_STRIDE * sizeof (v_t))) {
		MUL2_SETUP();

		XOR(REC_Q_X, REC_Q_X);
		for (c = firstdc; c < nbigcols; c++) {
			MUL2(REC_Q_X);
			if (c != x) {
				col = &rm->rm_col[c];
				XOR_ACC(COL_OFF(col, ioff), REC_Q_X);
			}
		}
		for (; c < ncols; c++)
			MUL2(REC_Q_X);

		XOR_ACC(COL_OFF(qcol, ioff), REC_Q_X);
		MUL(coeff[MUL_Q_X], REC_Q_X);
		STORE(COL_OFF(xcol, ioff), REC_Q_X);
	}
}

/*
 * Reconstruct single data column using Q parity
 * @rec_method	REC_Q_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */
static raidz_inline int
raidz_reconstruct_q_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;
	unsigned coeff[MUL_CNT];

	raidz_init_rec_coeff(rm, tgtidx, RAIDZ_REC_Q, coeff);

	raidz_math_begin();

	/* lastcsize */
	REC_Q_BLOCK(rm, 0, lastcsize, tgtidx, coeff, ncols, ncols);

	/* fullcols */
	if (lastcsize < xsize)
		REC_Q_BLOCK(rm, lastcsize, xsize, tgtidx, coeff,
		    ncols, nbigcols);

	raidz_math_end();

	return (1 << CODE_Q);
}


static raidz_inline void
REC_R_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const unsigned *coeff, const size_t ncols,
    const size_t nbigcols)
{
	int c;
	size_t ioff = 0;
	const int x = tgtidx[TARGET_X];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const rcol = &rm->rm_col[CODE_R];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t *col;

	REC_R_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_R_STRIDE * sizeof (v_t))) {
		MUL2_SETUP();
		XOR(REC_R_X, REC_R_X);
		for (c = firstdc; c < nbigcols; c++) {
			MUL4(REC_R_X);
			if (c != x) {
				col = &rm->rm_col[c];
				XOR_ACC(COL_OFF(col, ioff), REC_R_X);
			}
		}
		for (; c < ncols; c++)
			MUL4(REC_R_X);

		XOR_ACC(COL_OFF(rcol, ioff), REC_R_X);
		MUL(coeff[MUL_R_X], REC_R_X);
		STORE(COL_OFF(xcol, ioff), REC_R_X);
	}
}

/*
 * Reconstruct single data column using R parity
 * @rec_method	REC_R_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */
static raidz_inline int
raidz_reconstruct_r_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;
	unsigned coeff[MUL_CNT];

	raidz_init_rec_coeff(rm, tgtidx, RAIDZ_REC_R, coeff);

	raidz_math_begin();

	/* lastcsize */
	REC_R_BLOCK(rm, 0, lastcsize, tgtidx, coeff, ncols, ncols);

	/* fullcols */
	if (lastcsize < xsize)
		REC_R_BLOCK(rm, lastcsize, xsize, tgtidx, coeff,
		    ncols, nbigcols);

	raidz_math_end();

	return (1 << CODE_R);
}


static raidz_inline void
REC_PQ_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const unsigned *coeff, const size_t ncols,
    const size_t nbigcols, const boolean_t calcy)
{
	int c;
	size_t ioff = 0;
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const pcol = &rm->rm_col[CODE_P];
	raidz_col_t * const qcol = &rm->rm_col[CODE_Q];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t * const ycol = &rm->rm_col[y];
	raidz_col_t *col;

	REC_PQ_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_PQ_STRIDE * sizeof (v_t))) {
		LOAD(COL_OFF(pcol, ioff), REC_PQ_X);
		XOR(REC_PQ_Y, REC_PQ_Y);
		MUL2_SETUP();
		for (c = firstdc; c < nbigcols; c++) {
			MUL2(REC_PQ_Y);
			if (c != x && c != y) {
				col = &rm->rm_col[c];
				LOAD(COL_OFF(col, ioff), REC_PQ_D);
				XOR(REC_PQ_D, REC_PQ_X);
				XOR(REC_PQ_D, REC_PQ_Y);
			}
		}
		for (; c < ncols; c++)
			MUL2(REC_PQ_Y);

		XOR_ACC(COL_OFF(qcol, ioff), REC_PQ_Y);

		/* Save Pxy */
		COPY(REC_PQ_X, REC_PQ_D);

		/* Calc X */
		MUL(coeff[MUL_PQ_X], REC_PQ_X);
		MUL(coeff[MUL_PQ_Y], REC_PQ_Y);
		XOR(REC_PQ_Y,  REC_PQ_X);
		STORE(COL_OFF(xcol, ioff), REC_PQ_X);

		if (calcy) {
			/* Calc Y */
			XOR(REC_PQ_D,  REC_PQ_X);
			STORE(COL_OFF(ycol, ioff), REC_PQ_X);
		}
	}
}

/*
 * Reconstruct two data columns using PQ parity
 * @rec_method	REC_PQ_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */
static raidz_inline int
raidz_reconstruct_pq_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t ysize = rm->rm_col[y].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;
	unsigned coeff[MUL_CNT];

	raidz_init_rec_coeff(rm, tgtidx, RAIDZ_REC_PQ, coeff);

	raidz_math_begin();

	/* lastcsize */
	REC_PQ_BLOCK(rm, 0, lastcsize, tgtidx, coeff, ncols, ncols, B_TRUE);

	/* fullcols */
	if (lastcsize < ysize)
		REC_PQ_BLOCK(rm, lastcsize, ysize, tgtidx, coeff,
		    ncols, nbigcols, B_TRUE);

	if (ysize < xsize)
		REC_PQ_BLOCK(rm, ysize, xsize, tgtidx, coeff,
		    ncols, nbigcols, B_FALSE);

	raidz_math_end();

	return ((1 << CODE_P) | (1 << CODE_Q));
}


static raidz_inline void
REC_PR_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const unsigned *coeff, const size_t ncols,
    const size_t nbigcols, const boolean_t calcy)
{
	int c;
	size_t ioff;
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const pcol = &rm->rm_col[CODE_P];
	raidz_col_t * const rcol = &rm->rm_col[CODE_R];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t * const ycol = &rm->rm_col[y];
	raidz_col_t *col;

	REC_PR_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_PR_STRIDE * sizeof (v_t))) {
		LOAD(COL_OFF(pcol, ioff), REC_PR_X);
		XOR(REC_PR_Y, REC_PR_Y);
		MUL2_SETUP();
		for (c = firstdc; c < nbigcols; c++) {
			MUL4(REC_PR_Y);
			if (c != x && c != y) {
				col = &rm->rm_col[c];
				LOAD(COL_OFF(col, ioff), REC_PR_D);
				XOR(REC_PR_D, REC_PR_X);
				XOR(REC_PR_D, REC_PR_Y);
			}
		}
		for (; c < ncols; c++)
			MUL4(REC_PR_Y);

		XOR_ACC(COL_OFF(rcol, ioff), REC_PR_Y);

		/* Save Pxy */
		COPY(REC_PR_X,  REC_PR_D);

		/* Calc X */
		MUL(coeff[MUL_PR_X], REC_PR_X);
		MUL(coeff[MUL_PR_Y], REC_PR_Y);
		XOR(REC_PR_Y,  REC_PR_X);
		STORE(COL_OFF(xcol, ioff), REC_PR_X);

		if (calcy) {
			/* Calc Y */
			XOR(REC_PR_D,  REC_PR_X);
			STORE(COL_OFF(ycol, ioff), REC_PR_X);
		}
	}
}


/*
 * Reconstruct two data columns using PR parity
 * @rec_method	REC_PR_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */
static raidz_inline int
raidz_reconstruct_pr_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t ysize = rm->rm_col[y].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;
	unsigned coeff[MUL_CNT];

	raidz_init_rec_coeff(rm, tgtidx, RAIDZ_REC_PR, coeff);

	raidz_math_begin();

	/* lastcsize */
	REC_PR_BLOCK(rm, 0, lastcsize, tgtidx, coeff, ncols, ncols, B_TRUE);

	/* fullcols */
	if (lastcsize < ysize)
		REC_PR_BLOCK(rm, lastcsize, ysize, tgtidx, coeff,
		    ncols, nbigcols, B_TRUE);

	if (ysize < xsize)
		REC_PR_BLOCK(rm, ysize, xsize, tgtidx, coeff,
		    ncols, nbigcols, B_FALSE);

	raidz_math_end();

	return ((1 << CODE_P) | (1 << CODE_R));
}

static raidz_inline void
REC_QR_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const unsigned *coeff, const size_t ncols,
    const size_t nbigcols, const boolean_t calcy)
{
	int c;
	size_t ioff;
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const qcol = &rm->rm_col[CODE_Q];
	raidz_col_t * const rcol = &rm->rm_col[CODE_R];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t * const ycol = &rm->rm_col[y];
	raidz_col_t *col;

	REC_QR_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_QR_STRIDE * sizeof (v_t))) {
		MUL2_SETUP();
		XOR(REC_QR_X, REC_QR_X);
		XOR(REC_QR_Y, REC_QR_Y);
		for (c = firstdc; c < nbigcols; c++) {
			MUL2(REC_QR_X);
			MUL4(REC_QR_Y);
			if (c != x && c != y) {
				col = &rm->rm_col[c];
				LOAD(COL_OFF(col, ioff), REC_QR_D);
				XOR(REC_QR_D, REC_QR_X);
				XOR(REC_QR_D, REC_QR_Y);
			}
		}
		for (; c < ncols; c++) {
			MUL2(REC_QR_X);
			MUL4(REC_QR_Y);
		}

		XOR_ACC(COL_OFF(qcol, ioff), REC_QR_X);
		XOR_ACC(COL_OFF(rcol, ioff), REC_QR_Y);

		/* Save Qxy */
		COPY(REC_QR_X,  REC_QR_D);

		/* Calc X */
		MUL(coeff[MUL_QR_XQ], REC_QR_X);	/* X = Q * xqm */
		XOR(REC_QR_Y, REC_QR_X);		/* X = R ^ X   */
		MUL(coeff[MUL_QR_X], REC_QR_X);		/* X = X * xm  */
		STORE(COL_OFF(xcol, ioff), REC_QR_X);

		if (calcy) {
			/* Calc Y */
			MUL(coeff[MUL_QR_YQ], REC_QR_D); /* X = Q * xqm */
			XOR(REC_QR_Y, REC_QR_D);	 /* X = R ^ X   */
			MUL(coeff[MUL_QR_Y], REC_QR_D);	 /* X = X * xm  */
			STORE(COL_OFF(ycol, ioff), REC_QR_D);
		}
	}
}

/*
 * Reconstruct two data columns using QR parity
 * @rec_method	REC_QR_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */
static raidz_inline int
raidz_reconstruct_qr_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t ysize = rm->rm_col[y].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;
	unsigned coeff[MUL_CNT];

	raidz_init_rec_coeff(rm, tgtidx, RAIDZ_REC_QR, coeff);

	raidz_math_begin();

	/* lastcsize */
	REC_QR_BLOCK(rm, 0, lastcsize, tgtidx, coeff, ncols, ncols, B_TRUE);

	/* fullcols */
	if (lastcsize < ysize)
		REC_QR_BLOCK(rm, lastcsize, ysize, tgtidx, coeff,
		    ncols, nbigcols, B_TRUE);

	if (ysize < xsize)
		REC_QR_BLOCK(rm, ysize, xsize, tgtidx, coeff,
		    ncols, nbigcols, B_FALSE);

	raidz_math_end();

	return ((1 << CODE_Q) | (1 << CODE_R));
}


static raidz_inline void
REC_PQR_BLOCK(raidz_map_t * const rm, const size_t off, const size_t end,
    const int *tgtidx, const unsigned *coeff, const size_t ncols,
    const size_t nbigcols, const boolean_t calcy, const boolean_t calcz)
{
	int c;
	size_t ioff;
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const int z = tgtidx[TARGET_Z];
	const size_t firstdc = raidz_parity(rm);
	raidz_col_t * const pcol = &rm->rm_col[CODE_P];
	raidz_col_t * const qcol = &rm->rm_col[CODE_Q];
	raidz_col_t * const rcol = &rm->rm_col[CODE_R];
	raidz_col_t * const xcol = &rm->rm_col[x];
	raidz_col_t * const ycol = &rm->rm_col[y];
	raidz_col_t * const zcol = &rm->rm_col[z];
	raidz_col_t *col;

	REC_PQR_DEFINE();

	for (ioff = off; ioff < end; ioff += (REC_PQR_STRIDE * sizeof (v_t))) {
		MUL2_SETUP();
		LOAD(COL_OFF(pcol, ioff), REC_PQR_X);
		XOR(REC_PQR_Y, REC_PQR_Y);
		XOR(REC_PQR_Z, REC_PQR_Z);
		for (c = firstdc; c < nbigcols; c++) {
			MUL2(REC_PQR_Y);
			MUL4(REC_PQR_Z);
			if (c != x && c != y && c != z) {
				col = &rm->rm_col[c];
				LOAD(COL_OFF(col, ioff), REC_PQR_D);
				XOR(REC_PQR_D, REC_PQR_X);
				XOR(REC_PQR_D, REC_PQR_Y);
				XOR(REC_PQR_D, REC_PQR_Z);
			}
		}
		for (; c < ncols; c++) {
			MUL2(REC_PQR_Y);
			MUL4(REC_PQR_Z);
		}

		XOR_ACC(COL_OFF(qcol, ioff), REC_PQR_Y);
		XOR_ACC(COL_OFF(rcol, ioff), REC_PQR_Z);

		/* Save Pxyz and Qxyz */
		COPY(REC_PQR_X, REC_PQR_XS);
		COPY(REC_PQR_Y, REC_PQR_YS);

		/* Calc X */
		MUL(coeff[MUL_PQR_XP], REC_PQR_X);	/* Xp = Pxyz * xp   */
		MUL(coeff[MUL_PQR_XQ], REC_PQR_Y);	/* Xq = Qxyz * xq   */
		XOR(REC_PQR_Y, REC_PQR_X);
		MUL(coeff[MUL_PQR_XR], REC_PQR_Z);	/* Xr = Rxyz * xr   */
		XOR(REC_PQR_Z, REC_PQR_X);		/* X = Xp + Xq + Xr */
		STORE(COL_OFF(xcol, ioff), REC_PQR_X);

		if (calcy) {
			/* Calc Y */
			XOR(REC_PQR_X, REC_PQR_XS);	   /* Pyz = Pxyz + X */
			MUL(coeff[MUL_PQR_YU], REC_PQR_X); /* Xq = X * upd_q */
			XOR(REC_PQR_X, REC_PQR_YS);	   /* Qyz = Qxyz + Xq */
			COPY(REC_PQR_XS, REC_PQR_X);	   /* restore Pyz */
			MUL(coeff[MUL_PQR_YP], REC_PQR_X); /* Yp = Pyz * yp */
			MUL(coeff[MUL_PQR_YQ], REC_PQR_YS); /* Yq = Qyz * yq */
			XOR(REC_PQR_X, REC_PQR_YS);	    /* Y = Yp + Yq */
			STORE(COL_OFF(ycol, ioff), REC_PQR_YS);
		}

		if (calcy && calcz) {
			/* Calc Z */
			XOR(REC_PQR_XS, REC_PQR_YS);	/* Z = Pz = Pyz + Y */
			STORE(COL_OFF(zcol, ioff), REC_PQR_YS);
		}
	}
}

/*
 * Reconstruct three data columns using PQR parity
 * @rec_method	REC_PQR_BLOCK()
 *
 * @rm		RAIDZ map
 * @tgtidx	array of missing data indexes
 */
static raidz_inline int
raidz_reconstruct_pqr_impl(raidz_map_t *rm, const int *tgtidx)
{
	const int x = tgtidx[TARGET_X];
	const int y = tgtidx[TARGET_Y];
	const int z = tgtidx[TARGET_Z];
	const size_t ncols = raidz_ncols(rm);
	const size_t nbigcols = raidz_nbigcols(rm);
	const size_t xsize = rm->rm_col[x].rc_size;
	const size_t ysize = rm->rm_col[y].rc_size;
	const size_t zsize = rm->rm_col[z].rc_size;
	const size_t lastcsize = rm->rm_col[ncols-1].rc_size;
	unsigned coeff[MUL_CNT];

	raidz_init_rec_coeff(rm, tgtidx, RAIDZ_REC_PQR, coeff);

	raidz_math_begin();

	/* lastcsize */
	REC_PQR_BLOCK(rm, 0, lastcsize, tgtidx, coeff,
	    ncols, ncols, B_TRUE, B_TRUE);

	/* fullcols */
	if (lastcsize < zsize)
		REC_PQR_BLOCK(rm, lastcsize, zsize, tgtidx, coeff,
		    ncols, nbigcols, B_TRUE, B_TRUE);

	if (zsize < ysize)
		REC_PQR_BLOCK(rm, zsize, ysize, tgtidx, coeff,
		    ncols, nbigcols, B_TRUE, B_FALSE);

	if (ysize < xsize)
		REC_PQR_BLOCK(rm, ysize, xsize, tgtidx, coeff,
		    ncols, nbigcols, B_FALSE, B_FALSE);

	raidz_math_end();

	return ((1 << CODE_P) | (1 << CODE_Q) | (1 << CODE_R));
}

#endif /* _VDEV_RAIDZ_MATH_IMPL_H */
