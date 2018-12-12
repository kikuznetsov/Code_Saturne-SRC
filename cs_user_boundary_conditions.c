/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* Code_Saturne version 5.2.0 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_elec_model.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/
#include "cs_turbulence_model.h"
#include "cs_prototypes.h"

#include "read_from_rije_profile.h"
#include "read_from_ke_profile.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

#define Z0CABLE 0.001
#define Z0SEABED 0.0001
#define NUMOFLINES 120 //!!!! be carefull potential SIGSEV!!! should be less than numbers of lines in the file
#define FILEPROFILE "tmpUx.csv"
//#define FILEPROFILE "/home/konst/Projects/STHYF/Calcs/current-cylinder-bc/INIT/Ux.csv"
/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in]     nvar          total number of variable BC's
 * \param[in]     bc_type       boundary face types
 * \param[in]     icodcl        boundary face code
 *                                - 1  -> Dirichlet
 *                                - 2  -> convective outlet
 *                                - 3  -> flux density
 *                                - 4  -> sliding wall and u.n=0 (velocity)
 *                                - 5  -> friction and u.n=0 (velocity)
 *                                - 6  -> roughness and u.n=0 (velocity)
 *                                - 9  -> free inlet/outlet (velocity)
 *                                inflowing possibly blocked
 * \param[in]     rcodcl        boundary condition values
 *                                rcodcl(3) = flux density value
 *                                (negative for gain) in W/m2
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(int         nvar,
                            int         bc_type[],
                            int         icodcl[],
                            cs_real_t   rcodcl[])
{
  //PREPARE NON-INLET BOUNDARIES
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *) cs_glob_mesh_quantities->b_face_normal;
  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;
  cs_lnum_t *lstelt = NULL;
  cs_lnum_t *face_list;
  cs_lnum_t  nelts;
  cs_lnum_t  cell_id;
  cs_field_t *f;
  cs_field_t *fu = CS_F_(u);
  const int keyvar = cs_field_key_id("variable_id");

  int ivar_Ux = cs_field_get_key_int(fu, keyvar) - 1 + 0;  //var for Ux
  int ivar_Uy = cs_field_get_key_int(fu, keyvar) - 1 + 1;  //var for Uy
  int ivar_Uz = cs_field_get_key_int(fu, keyvar) - 1 + 2;  //var for Uz
  //const int keyRough = cs_field_key_id("boundary_roughness");

  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);
  
  //top
  cs_selector_get_b_face_list("top", &nelts, lstelt);
  for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
    cs_lnum_t face_id = lstelt[ilelt];
    bc_type[face_id] = CS_SYMMETRY;

  }

  //bottom
  // cs_selector_get_b_face_list("bottom", &nelts, lstelt);
  // for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
  //   cs_lnum_t face_id = lstelt[ilelt];
  //   bc_type[face_id] = CS_SYMMETRY;
  // }
  // cs_selector_get_b_face_list("bottom", &nelts, lstelt);
  // //int ivarv = cs_field_get_key_int(fu, keyvar) - 1; // access to U_x 
  // for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
  //   cs_lnum_t face_id = lstelt[ilelt];
  //   bc_type[face_id] = CS_ROUGHWALL;
  //   rcodcl[2*n_b_faces*nvar + ivar_Ux*n_b_faces + face_id] = Z0SEABED; //scalar roguhness
  //   //rcodcl[2*n_b_faces*nvar + ivar_Uy*n_b_faces + face_id] = Z0SEABED; //roughness value z0
  // }

  //for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
//    cs_lnum_t face_id = lstelt[ilelt];
//    bc_type[face_id] = CS_OUTLET;
//    cell_id = b_face_cells[face_id]; // associated boundary cell
//}

  ///////////PREPARE INLET BOUNDARIES

  ////////prepare input files
  //define file name of profile and length
  size_t num_lines = NUMOFLINES; //number of points in profile defined by user in GUI
  const char* fName = FILEPROFILE;

  ///IF k-epsilon models
  if( cs_glob_turb_model->itytur==2){///k-epsilon
    //printf("k-epsilon\n");

    cs_real_t *k = (cs_real_t *)(CS_F_(k)->val); 

    //Prepare array for inital profiles
    struct profile_keps_t* profile =
                      (struct profile_keps_t *) malloc(sizeof(struct profile_keps_t)
                      + num_lines*sizeof(struct record_keps_t));
    profile->n_rows = num_lines;
    profile->rec = (struct record_keps_t *) malloc ((profile->n_rows)*sizeof(struct record_keps_t));

    //Read csv file
    int status = read_profile_keps(fName, num_lines, profile);
    if(status==EXIT_FAILURE){
      printf("error of reading file\n");
      return;
    }

    // Inlet
    cs_real_t y_curr;
    struct record_keps_t temp;
    f = CS_F_(k);//get turbulence energy field
    int ivar_k = cs_field_get_key_int(f, keyvar) - 1;  //ivar for k
    
    f = CS_F_(eps);//get turbulence energy field
    int ivar_eps = cs_field_get_key_int(f, keyvar) - 1;  //ivar for eps


    cs_selector_get_b_face_list("inlet or outlet", &nelts, lstelt);
    for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
      cs_lnum_t face_id = lstelt[ilelt];
      bc_type[face_id] = CS_INLET;
      cell_id = b_face_cells[face_id]; // associated boundary cell
      
      y_curr = cell_cen[cell_id][1]; // define height
      temp = interpolate_keps(profile, y_curr);

      //define Ux
      icodcl[(ivar_Ux) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_Ux) * n_b_faces + face_id] = temp.u; //Value 

      //define Uy
      icodcl[(ivar_Uy) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_Uy) * n_b_faces + face_id] = temp.v; //Value  

      //define Uz
      icodcl[(ivar_Uz) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_Uz) * n_b_faces + face_id] = 0.; //Value  

      //define k  
      icodcl[(ivar_k) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_k) * n_b_faces + face_id] = temp.k; //Value  

      //define eps
      icodcl[(ivar_eps) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_eps) * n_b_faces + face_id] = temp.eps; //Value  

    }
    //Deallocate the memory
    free(profile->rec);
    free(profile);

  }
  ///IF Rij-epsilon models (SSG,LRR,EBRSM)
  else if( cs_glob_turb_model->itytur==3){//Rij-epsilon

    //cs_real_6_t *rij = (cs_real_6_t *)(CS_F_(rij)->val); 

    printf("SSG\n"); 
    //Prepare array for inital profiles
    struct profile_rijssg_t* profile =
                      (struct profile_rijssg_t *) malloc(sizeof(struct profile_rijssg_t)
                      + num_lines*sizeof(struct record_rijssg_t));;
    profile->n_rows = num_lines;
    profile->rec = (struct record_rijssg_t *) malloc ((profile->n_rows)*sizeof(struct record_rijssg_t));

    //Read csv file
    int status = read_profile_SSG(fName, num_lines, profile);
    if(status==EXIT_FAILURE){
      printf("error of reading file\n");
      return;
    }

    // Inlet
    cs_real_t y_curr;
    struct record_rijssg_t temp;
    f = CS_F_(rij);//get turbulence energy field
    int ivar_R0 = cs_field_get_key_int(f, keyvar) - 1 + 0;  //var for R_xx
  	int ivar_R1 = cs_field_get_key_int(f, keyvar) - 1 + 1;  //var for R_yy
  	int ivar_R2 = cs_field_get_key_int(f, keyvar) - 1 + 2;  //var for R_zz
  	int ivar_R3 = cs_field_get_key_int(f, keyvar) - 1 + 3;  //var for R_xy
  	int ivar_R4 = cs_field_get_key_int(f, keyvar) - 1 + 4;  //var for R_yz
  	int ivar_R5 = cs_field_get_key_int(f, keyvar) - 1 + 5;  //var for R_xz
    
    f = CS_F_(eps);//get turbulence energy field
    int ivar_eps = cs_field_get_key_int(f, keyvar) - 1;  //ivar for eps

    //Inlet or Outlet
    cs_selector_get_b_face_list("inlet or outlet", &nelts, lstelt);
    for (cs_lnum_t ilelt = 0; ilelt < nelts; ilelt++) {
      cs_lnum_t face_id = lstelt[ilelt];
      bc_type[face_id] = CS_INLET;
      cell_id = b_face_cells[face_id]; // associated boundary cell
      
      y_curr = cell_cen[cell_id][1]; // define height
      temp = interpolate_rijssg(profile, y_curr);

      //define Ux
      icodcl[(ivar_Ux) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_Ux) * n_b_faces + face_id] = temp.u; //Value 

      //define Uy
      icodcl[(ivar_Uy) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_Uy) * n_b_faces + face_id] = temp.v; //Value  

      //define Uz
      icodcl[(ivar_Uz) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_Uz) * n_b_faces + face_id] = 0.; //Value  

      //define R_xx  
      icodcl[(ivar_R0) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_R0) * n_b_faces + face_id] = temp.rxx; //Value  
      //define R_yy  
      icodcl[(ivar_R1) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_R1) * n_b_faces + face_id] = temp.ryy; //Value  
      //define R_zz  
      icodcl[(ivar_R2) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_R2) * n_b_faces + face_id] = temp.rzz; //Value  
      //define R_xy
      icodcl[(ivar_R3) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_R3) * n_b_faces + face_id] = temp.rxy; //Value  
      //define R_yz  
      icodcl[(ivar_R4) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_R4) * n_b_faces + face_id] = temp.ryz; //Value  
      //define R_xz  
      icodcl[(ivar_R5) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_R5) * n_b_faces + face_id] = temp.rxz; //Value  

      //define eps
      icodcl[(ivar_eps) * n_b_faces + face_id] = 1; //Dirihlet value
      rcodcl[(ivar_eps) * n_b_faces + face_id] = temp.eps; //Value  

    }
    //Deallocate the memory
    free(profile->rec);
    free(profile);
  }


  BFT_FREE(lstelt);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
