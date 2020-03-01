/*    
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.   
  Do not edit this file directly as it can be overwritten!!   
   
  This Source Code Form is subject to the terms of the Mozilla Public   
  License, v. 2.0. If a copy of the MPL was not distributed with this   
  file, You can obtain one at http://mozilla.org/MPL/2.0/.   
   
  Maple version     : Maple 2016 (X86 64 LINUX)   
  Maple source      : ../maple/gga_c_w94.mpl   
  Type of functional: work_gga_c   
*/   
   
#ifdef DEVICE   
__device__ void xc_gga_c_w94_func   
  (const void *p, xc_gga_work_c_t *r)   
#else   
void xc_gga_c_w94_func   
  (const xc_func_type *p, xc_gga_work_c_t *r)   
#endif   
{   
  double t1, t4, t5, t6, t8, t9, t10, t12;   
  double t13, t14, t17, t18, t19, t20, t22, t23;   
  double t24, t28, t29, t31, t32, t33, t36, t37;   
  double t38, t40, t43, t44, t48, t52, t54, t55;   
  double t56, t59, t60, t68, t73, t74, t75, t76;   
  double t77, t80, t81, t85, t88, t89, t92, t95;   
  double t103, t105, t107, t110, t114, t154, t158, t159;   
  double t179;   
   
   
  t1 = Heaviside(-r->z);   
  t4 = -0.2e1 * r->z * t1 + r->z;   
  t5 = cbrt(t4);   
  t6 = t5 * t5;   
  t8 = -t6 * t4 + 0.1e1;   
  t9 = sqrt(t8);   
  t10 = r->xt * r->xt;   
  t12 = pow(r->xt, 0.1e1 / 0.16e2);   
  t13 = t12 * t12;   
  t14 = t13 * t12;   
  t17 = M_CBRT3;   
  t18 = t17 * t17;   
  t19 = M_CBRT4;   
  t20 = t18 * t19;   
  t22 = cbrt(0.1e1 / 0.31415926535897932385e1);   
  t23 = 0.1e1 / t22;   
  t24 = t23 * t10;   
  t28 = 0.118e2 + 0.150670e0 * t14 * t10 * r->xt + 0.36733333333333333333e-2 * t20 * t24 * r->rs + r->rs;   
  t29 = 0.1e1 / t28;   
  r->f = -t9 * t29;   
   
  if(r->order < 1) return;   
   
  t31 = t28 * t28;   
  t32 = 0.1e1 / t31;   
  t33 = t9 * t32;   
  t36 = 0.36733333333333333333e-2 * t20 * t24 + 0.1e1;   
  r->dfdrs = t33 * t36;   
  t37 = 0.1e1 / t9;   
  t38 = t37 * t29;   
  t40 = 0.0;   
  t43 = 0.2e1 * r->z * t40 - 0.2e1 * t1 + 0.1e1;   
  t44 = t6 * t43;   
  r->dfdz = 0.5e1 / 0.6e1 * t38 * t44;   
  t48 = t23 * r->xt;   
  t52 = 0.48026062500000000000e0 * t14 * t10 + 0.73466666666666666666e-2 * t20 * t48 * r->rs;   
  r->dfdxt = t33 * t52;   
  r->dfdxs[0] = 0.0e0;   
  r->dfdxs[1] = 0.0e0;   
   
  if(r->order < 2) return;   
   
  t54 = 0.1e1 / t31 / t28;   
  t55 = t9 * t54;   
  t56 = t36 * t36;   
  r->d2fdrs2 = -0.2e1 * t55 * t56;   
  t59 = t37 * t32;   
  t60 = t36 * t6;   
  r->d2fdrsz = -0.5e1 / 0.6e1 * t59 * t60 * t43;   
  t68 = t19 * t23;   
  r->d2fdrsxt = -0.2e1 * t55 * t36 * t52 + 0.73466666666666666666e-2 * t33 * t18 * t68 * r->xt;   
  r->d2fdrsxs[0] = 0.0e0;   
  r->d2fdrsxs[1] = 0.0e0;   
  t73 = 0.1e1 / t9 / t8;   
  t74 = t73 * t29;   
  t75 = t5 * t4;   
  t76 = t43 * t43;   
  t77 = t75 * t76;   
  t80 = 0.1e1 / t5;   
  t81 = t80 * t76;   
  t85 = 0.0;   
  t88 = 0.2e1 * r->z * t85 + 0.4e1 * t40;   
  t89 = t6 * t88;   
  r->d2fdz2 = 0.25e2 / 0.36e2 * t74 * t77 + 0.5e1 / 0.9e1 * t38 * t81 + 0.5e1 / 0.6e1 * t38 * t89;   
  t92 = t44 * t52;   
  r->d2fdzxt = -0.5e1 / 0.6e1 * t59 * t92;   
  r->d2fdzxs[0] = 0.0e0;   
  r->d2fdzxs[1] = 0.0e0;   
  t95 = t52 * t52;   
  t103 = 0.10505701171875000000e1 * t14 * r->xt + 0.73466666666666666666e-2 * t20 * t23 * r->rs;   
  r->d2fdxt2 = t33 * t103 - 0.2e1 * t55 * t95;   
  r->d2fdxtxs[0] = 0.0e0;   
  r->d2fdxtxs[1] = 0.0e0;   
  r->d2fdxs2[0] = 0.0e0;   
  r->d2fdxs2[1] = 0.0e0;   
  r->d2fdxs2[2] = 0.0e0;   
   
  if(r->order < 3) return;   
   
  t105 = t31 * t31;   
  t107 = t9 / t105;   
  r->d3fdrs3 = 0.6e1 * t107 * t56 * t36;   
  t110 = t37 * t54;   
  r->d3fdrs2z = 0.5e1 / 0.3e1 * t110 * t56 * t6 * t43;   
  t114 = t73 * t32;   
  r->d3fdrsz2 = -0.25e2 / 0.36e2 * t114 * t36 * t75 * t76 - 0.5e1 / 0.9e1 * t59 * t36 * t80 * t76 - 0.5e1 / 0.6e1 * t59 * t60 * t88;   
  r->d3fdrszxt = 0.5e1 / 0.3e1 * t110 * t36 * t92 - 0.61222222222222222222e-2 * t59 * t20 * t48 * t44;   
  r->d3fdrszxs[0] = 0.0e0;   
  r->d3fdrszxs[1] = 0.0e0;   
  r->d3fdrs2xt = 0.6e1 * t107 * t56 * t52 - 0.29386666666666666666e-1 * t55 * t36 * t20 * t48;   
  r->d3fdrsxt2 = 0.6e1 * t107 * t36 * t95 - 0.29386666666666666666e-1 * t55 * t18 * t68 * r->xt * t52 - 0.2e1 * t55 * t36 * t103 + 0.73466666666666666666e-2 * t33 * t20 * t23;   
  r->d3fdrsxtxs[0] = 0.0e0;   
  r->d3fdrsxtxs[1] = 0.0e0;   
  r->d3fdrs2xs[0] = 0.0e0;   
  r->d3fdrs2xs[1] = 0.0e0;   
  r->d3fdrsxs2[0] = 0.0e0;   
  r->d3fdrsxs2[1] = 0.0e0;   
  r->d3fdrsxs2[2] = 0.0e0;   
  t154 = t8 * t8;   
  t158 = t4 * t4;   
  t159 = t76 * t43;   
  t179 = 0.0;   
  r->d3fdz3 = 0.125e3 / 0.72e2 / t9 / t154 * t29 * t158 * t159 + 0.25e2 / 0.18e2 * t74 * t5 * t159 + 0.25e2 / 0.12e2 * t74 * t75 * t43 * t88 - 0.5e1 / 0.27e2 * t38 / t75 * t159 + 0.5e1 / 0.3e1 * t38 * t80 * t43 * t88 + 0.5e1 / 0.6e1 * t38 * t6 * (0.2e1 * r->z * t179 + 0.6e1 * t85);   
  r->d3fdz2xt = -0.25e2 / 0.36e2 * t114 * t77 * t52 - 0.5e1 / 0.9e1 * t59 * t81 * t52 - 0.5e1 / 0.6e1 * t59 * t89 * t52;   
  r->d3fdzxt2 = 0.5e1 / 0.3e1 * t110 * t44 * t95 - 0.5e1 / 0.6e1 * t59 * t44 * t103;   
  r->d3fdzxtxs[0] = 0.0e0;   
  r->d3fdzxtxs[1] = 0.0e0;   
  r->d3fdz2xs[0] = 0.0e0;   
  r->d3fdz2xs[1] = 0.0e0;   
  r->d3fdzxs2[0] = 0.0e0;   
  r->d3fdzxs2[1] = 0.0e0;   
  r->d3fdzxs2[2] = 0.0e0;   
  r->d3fdxt3 = 0.6e1 * t107 * t95 * t52 - 0.6e1 * t55 * t52 * t103 + 0.12475520141601562500e1 * t33 * t14;   
  r->d3fdxt2xs[0] = 0.0e0;   
  r->d3fdxt2xs[1] = 0.0e0;   
  r->d3fdxtxs2[0] = 0.0e0;   
  r->d3fdxtxs2[1] = 0.0e0;   
  r->d3fdxtxs2[2] = 0.0e0;   
  r->d3fdxs3[0] = 0.0e0;   
  r->d3fdxs3[1] = 0.0e0;   
  r->d3fdxs3[2] = 0.0e0;   
  r->d3fdxs3[3] = 0.0e0;   
   
  if(r->order < 4) return;   
   
   
}   
   
#ifndef DEVICE   
#define maple2c_order 3   
#define maple2c_func  xc_gga_c_w94_func   
#define kernel_id 26 
#endif   