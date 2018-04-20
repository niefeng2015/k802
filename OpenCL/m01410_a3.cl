/**
 * Author......: See docs/credits.txt
 * License.....: MIT
 */

#define NEW_SIMD_CODE

#include "inc_vendor.cl"
#include "inc_hash_constants.h"
#include "inc_hash_functions.cl"
#include "inc_types.cl"
#include "inc_common.cl"
#include "inc_simd.cl"
#include "inc_hash_sha256.cl"


__kernel void m01410_mxx (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const pw_t *combs_buf, __constant const u32x *words_buf_r, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * modifier
   */

  const u32 lid = get_local_id (0);
  const u32 gid = get_global_id (0);

  if (gid >= gid_max) return;

  /**
   * base
   */

  const u32 pw_len = pws[gid].pw_len;

  const u32 pw_lenv = ceil ((float) pw_len / 4);

  u32x w[64] = { 0 };

  for (int idx = 0; idx < pw_lenv; idx++)
  {
    w[idx] = pws[gid].i[idx];

    barrier (CLK_GLOBAL_MEM_FENCE);
  }

  const u32 salt_len = salt_bufs[salt_pos].salt_len;

  const u32 salt_lenv = ceil ((float) salt_len / 4);

  u32x s[64] = { 0 };

  for (int idx = 0; idx < salt_lenv; idx++)
  {
    s[idx] = swap32_S (salt_bufs[salt_pos].salt_buf[idx]);

    barrier (CLK_GLOBAL_MEM_FENCE);
  }

  /**
   * loop
   */

  u32x w0l = w[0];

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    const u32x w0r = words_buf_r[il_pos / VECT_SIZE];

    const u32x w0 = w0l | w0r;

    w[0] = w0;

    sha256_ctx_vector_t ctx;

    sha256_init_vector (&ctx);

    sha256_update_vector (&ctx, w, pw_len);

    sha256_update_vector (&ctx, s, salt_len);

    sha256_final_vector (&ctx);

    const u32x r0 = ctx.h[DGST_R0];
    const u32x r1 = ctx.h[DGST_R1];
    const u32x r2 = ctx.h[DGST_R2];
    const u32x r3 = ctx.h[DGST_R3];

    COMPARE_M_SIMD (r0, r1, r2, r3);
  }
    //printf("hello000\n");
}

__kernel void m01410_sxx (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const pw_t *combs_buf, __constant const u32x *words_buf_r, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * modifier
   */

  const u32 lid = get_local_id (0);
  const u32 gid = get_global_id (0);

  if (gid >= gid_max) return;

  /**
   * digest
   */

  const u32 search[4] =
  {
    digests_buf[digests_offset].digest_buf[DGST_R0],
    digests_buf[digests_offset].digest_buf[DGST_R1],
    digests_buf[digests_offset].digest_buf[DGST_R2],
    digests_buf[digests_offset].digest_buf[DGST_R3]
  };

  //printf("[DEBUGINFO] search[4] = %08x %08x %08x %08x\n",search[0],search[1],search[2],search[3]);

  /**
   * base
   */

  const u32 pw_len = pws[gid].pw_len;

  const u32 pw_lenv = ceil ((float) pw_len / 4);

  u32x w[64] = { 0 };

  for (int idx = 0; idx < pw_lenv; idx++)
  {
    w[idx] = pws[gid].i[idx];

    barrier (CLK_GLOBAL_MEM_FENCE);
  }

  //printf("[DEBUGINFO] len=%d,pws = %08x %08x\n",pw_len,w[0],w[1]);

  const u32 salt_len = salt_bufs[salt_pos].salt_len;

  const u32 salt_lenv = ceil ((float) salt_len / 4);

  u32x s[64] = { 0 };

   //printf salt
   //printf("[DEBUGINFO]salt_len=%d,salt_lenv=%d,salt:[",salt_len,salt_lenv);
   //for (int idx = 0; idx < salt_lenv; idx++)
   //{
   //  printf("%08x ",salt_bufs[salt_pos].salt_buf[idx]);
   //}
   //printf("]\n");

  for (int idx = 0; idx < salt_lenv; idx++)
  {
    s[idx] = swap32_S (salt_bufs[salt_pos].salt_buf[idx]);
    barrier (CLK_GLOBAL_MEM_FENCE);
  }
//    printf("[DEBUGINFO]after swap32:\n%08x\n%08x\n",s[0],s[1]);


  /**
   * loop
   */

  u32x w0l = w[0];

  printf("[DEBUGINFO] il_cnt=%d,VECT_SIZE=%d\n",il_cnt,VECT_SIZE);
  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    const u32x w0r = words_buf_r[il_pos / VECT_SIZE];

    const u32x w0 = w0l | w0r;

    w[0] = w0;

    //printf("[DEBUGINFO]\n w[0]=%08x w[1]=%08x\n w0l=%08x, w0r=%08x, w0l|w0r=%08x\n",w0l,w[1],w0l,w0r,w0);
    printf("[w]\n%08x %08x\n [s]\n%08x %08x\n",w[0],w[1],s[0],s[1]);
	printf("pw_len=%d,salt_len=%d\n",pw_len,salt_len);

    sha256_ctx_vector_t ctx;

    sha256_init_vector (&ctx);

    sha256_update_vector (&ctx, w, pw_len);

    sha256_update_vector (&ctx, s, salt_len);

    sha256_final_vector (&ctx);

    const u32x r0 = ctx.h[DGST_R0];
    const u32x r1 = ctx.h[DGST_R1];
    const u32x r2 = ctx.h[DGST_R2];
    const u32x r3 = ctx.h[DGST_R3];

	printf("[DIGEST]\n%08x %08x %08x %08x\n[SEARCH]\n%08x %08x %08x %08x\n",r0,r1,r2,r3,search[0],search[1],search[2],search[3]);

    COMPARE_S_SIMD (r0, r1, r2, r3);
    //printf("sha256($pass.$salt) runs\n");
  }
}
