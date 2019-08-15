/* 
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ../maple/gga_x_2d_b88.mpl
  Type of functional: work_gga_x
*/

#ifdef DEVICE
__device__ void xc_gga_x_2d_b88_enhance
  (const void *p,  xc_gga_work_x_t *r)
#else
void xc_gga_x_2d_b88_enhance
  (const xc_func_type *p,  xc_gga_work_x_t *r)
#endif
{
  double t1, t2, t5, t6, t11, t12, t13, t15;
  double t16, t17, t20, t24, t28, t29, t30, t35;
  double t38, t48, t60;


  t1 = r->x * r->x;
  t2 = log(r->x + sqrt(r->x * r->x + 0.1e1));
  t5 = 0.1e1 + 0.149128e0 * r->x * t2;
  t6 = 0.1e1 / t5;
  r->f = 0.1e1 + 0.12390117088023646599e-1 * t1 * t6;

  if(r->order < 1) return;

  t11 = t5 * t5;
  t12 = 0.1e1 / t11;
  t13 = t1 * t12;
  t15 = t1 + 0.1e1;
  t16 = sqrt(t15);
  t17 = 0.1e1 / t16;
  t20 = 0.149128e0 * t2 + 0.149128e0 * r->x * t17;
  r->dfdx = 0.24780234176047293198e-1 * r->x * t6 - 0.12390117088023646599e-1 * t13 * t20;

  if(r->order < 2) return;

  t24 = r->x * t12;
  t28 = 0.1e1 / t11 / t5;
  t29 = t1 * t28;
  t30 = t20 * t20;
  t35 = 0.1e1 / t16 / t15;
  t38 = 0.298256e0 * t17 - 0.149128e0 * t1 * t35;
  r->d2fdx2 = 0.24780234176047293198e-1 * t6 - 0.49560468352094586396e-1 * t24 * t20 + 0.24780234176047293198e-1 * t29 * t30 - 0.12390117088023646599e-1 * t13 * t38;

  if(r->order < 3) return;

  t48 = t11 * t11;
  t60 = t15 * t15;
  r->d3fdx3 = -0.74340702528141879594e-1 * t12 * t20 + 0.14868140505628375919e0 * r->x * t28 * t30 - 0.74340702528141879594e-1 * t24 * t38 - 0.74340702528141879594e-1 * t1 / t48 * t30 * t20 + 0.74340702528141879594e-1 * t29 * t20 * t38 - 0.12390117088023646599e-1 * t13 * (-0.596512e0 * t35 * r->x + 0.447384e0 * t1 * r->x / t16 / t60);

  if(r->order < 4) return;


}

#ifndef DEVICE
#define maple2c_order 3
#define maple2c_func  xc_gga_x_2d_b88_enhance
#endif