/*============================================================================
 * User initialization prior to solving time steps.
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
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_lagr.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include "read_from_rije_profile.h"
#include "read_from_ke_profile.h"
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


#define NUMOFLINES_INIT 120
#define FILEPROFILE_INIT "tmpUx.csv"
//#define FILEPROFILE_INIT "/home/konst/Projects/STHYF/Calcs/current-cylinder-bc/INIT/Ux.csv"
/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(void)
{

  //define file name of profile and length
  //!!!! be carefull potential SIGSEV!!! should be less than numbers of lines in the file
  size_t num_lines = NUMOFLINES_INIT; //number of points in profile defined by user in GUI

  const char* fName = FILEPROFILE_INIT;

    // Define CS-variables
  const int location_id = CS_MESH_LOCATION_CELLS;
  const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0]; //Numbers of cells
  cs_real_3_t *vel = (cs_real_3_t *)(CS_F_(u)->val); //velocity field

  

  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  const cs_real_3_t* xyz = (const cs_real_3_t *)cs_glob_mesh->vtx_coord;

  cs_real_t *eps = (cs_real_t *)(CS_F_(eps)->val);
 
  
  ///IF k-epsilon models
  if( cs_glob_turb_model->itytur==2){///k-epsilon
    printf("k-epsilon\n");

    cs_real_t *k = (cs_real_t *)(CS_F_(k)->val); 

    //Prepare array for inital profiles
    struct profile_keps_t* profile =
                      (struct profile_keps_t *) malloc(sizeof(struct profile_keps_t)
                      + num_lines*sizeof(struct record_keps_t));;
    profile->n_rows = num_lines;
    profile->rec = (struct record_keps_t *) malloc ((profile->n_rows)*sizeof(struct record_keps_t));

    //Read csv file
    int status = read_profile_keps(fName, num_lines, profile);
    if(status==EXIT_FAILURE){
      printf("error of reading file\n");
      return;
    }
    cs_real_t y_curr;
    struct record_keps_t temp;
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      y_curr = cell_cen[i][1];
      temp = interpolate_keps(profile, y_curr);
      vel[i][0]=temp.u;
      vel[i][1]=temp.v;
      vel[i][2]=0.0;
      
      k[i] = temp.k;
      eps[i] = temp.eps;
    }
    //Deallocate the memory
    free(profile->rec);
    free(profile);

  }
  ///IF Rij-epsilon models (SSG,LRR,EBRSM)
  else if( cs_glob_turb_model->itytur==3){//Rij-epsilon

    cs_real_6_t *rij = (cs_real_6_t *)(CS_F_(rij)->val); 

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
    cs_real_t y_curr;
    struct record_rijssg_t temp;
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      y_curr = cell_cen[i][1];
      temp = interpolate_rijssg(profile, y_curr);
      vel[i][0]=temp.u;
      vel[i][1]=temp.v;
      vel[i][2]=0.0;
      
      eps[i] = temp.eps;

      rij[i][0] = temp.rxx; //R_xx
      rij[i][1] = temp.ryy; //R_yy
      rij[i][2] = temp.rzz; //R_zz

      rij[i][3] = temp.rxy;  //R_xy
      rij[i][4] = temp.ryz;  //R_yz
      rij[i][5] = temp.rxz;  //R_xz
    }
    //Deallocate the memory
    free(profile->rec);
    free(profile);
  }
  else{
    printf("Error!There is no user-defined initialization for that turbulence model!\n");
    return;
  }


}

/*----------------------------------------------------------------------------*/

END_C_DECLS
