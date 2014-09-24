////////////////////////////////////////////////////////////////////////////////
// create_healpix_from_cat.C
// Purpose: Create a Healpix map from a galaxy catalog
// Usage: create_healpix_from_cat [catalog_file in FITS format] [nside]
// Assumes catalog has RA, DEC columns
// Author: Nacho Sevilla
// Last modified: Sep 23 2014
// To compile (pcae75)
// gcc -o create_healpix_from_cat create_healpix_from_cat.C -L/usr/local/cfitsio/lib -lcfitsio -L/usr/local/healpix/lib -lchealpix -lstdc++ 
////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "/usr/local/cfitsio/include/fitsio.h"
#include "/usr/local/healpix/include/chealpix.h"

using namespace std;

int main(int argc, char *argv[]){

  float d2r = M_PI/180.0;
  int n;

  fitsfile *catalog_fptr;
  fitsfile *healpix_map_out_fptr;
  int sysmap_status=0;
  int catalog_status=0;
  int status=0;
  char expr[80];
  int catalog_hdunum,catalog_nhdu,catalog_ityphdu;
  int nmoverel=1,nmoveabs=2;
  long catalog_nrows;
  int catalog_ncols;
  int anynull;
  long longnull = 0;
  float fltnull = 0.0;
  
  int nside;
  long mapsize; 
  long *hp_map_value;
  long *pixel_nb;

  float *ra,*dec;
  char racolname[8]="ra";
  char deccolname[8]="dec";
  char pixcolname[8]="pixel";
  char ngalcolname[8]="ngal";
  int racolnum;
  int deccolnum;
  long pix;

  char tform_pix[8] = "D";
  char tform_ngal[8] = "D";

  if (argc != 3){
    printf("Usage: create_healpix_from_cat catalog_file nside\n");
    return 1;
  }

  nside = atoi(argv[2]);
  mapsize = 12*nside*nside;

  //allocate map
  hp_map_value = (long *)calloc(mapsize,sizeof(long));
  pixel_nb = (long *)calloc(mapsize,sizeof(long));

  //initialize map
  for(n=0;n<mapsize;n++){hp_map_value[n]=0;}
  for(n=0;n<mapsize;n++){pixel_nb[n]=n;}

  ////open catalog
  printf("\nCatalog: %s\n",argv[1]);
  if(!fits_open_table(&catalog_fptr,argv[1],READONLY,&catalog_status)){
    fits_get_num_hdus(catalog_fptr, &catalog_hdunum, &catalog_status);
    fits_get_hdu_num(catalog_fptr,&catalog_nhdu);
    fits_get_hdu_type(catalog_fptr,&catalog_ityphdu,&catalog_status);
    fits_movabs_hdu(catalog_fptr,nmoveabs,NULL,&catalog_status);
    fits_get_hdu_num(catalog_fptr,&catalog_nhdu);
    fits_get_num_rows(catalog_fptr,&catalog_nrows,&catalog_status);
    fits_get_num_cols(catalog_fptr,&catalog_ncols,&catalog_status);
    printf("Number of rows %ld\n",catalog_nrows);
    printf("Number of columns %ld\n",catalog_ncols);
  }
  printf("Catalog STATUS:%d\n",catalog_status);
  if(catalog_status!=0){
    exit(1);
  }
  ////read catalog
  // Assumes catalog has RA, DEC columns
  ra = (float *) malloc((int)catalog_nrows*sizeof(float));
  dec = (float *) malloc((int)catalog_nrows*sizeof(float));
  fits_get_colnum(catalog_fptr,CASEINSEN,racolname,&racolnum,&status);
  fits_get_colnum(catalog_fptr,CASEINSEN,deccolname,&deccolnum,&status);
  fits_read_col(catalog_fptr,TFLOAT,racolnum,1L,1L,catalog_nrows,&fltnull,ra,&anynull,&status);
  printf("Read column ra, status %d\n",status);
  fits_read_col(catalog_fptr,TFLOAT,deccolnum,1L,1L,catalog_nrows,&fltnull,dec,&anynull,&status);
  printf("Read column dec, status %d\n",status);
  
  ////for each row in the catalog, calculate theta and phi and write array with value corresponding
  //// to the pixel at that position
  for(n=0;n<catalog_nrows;n++){
    ang2pix_ring(nside,M_PI/2.0-dec[n]*d2r,ra[n]*d2r,&pix);
    hp_map_value[pix]++;
  }

  fits_create_file(&healpix_map_out_fptr,"!./out.fits",&status);
  printf("Created out.fits STATUS:%d\n",status);

  fits_create_tbl(healpix_map_out_fptr,BINARY_TBL,0,0,NULL,NULL,NULL,NULL,&status);
  printf("Created table STATUS:%d\n",status);  

  fits_insert_col(healpix_map_out_fptr,1,pixcolname,tform_pix,&status);
  printf("Inserted column %s STATUS:%d\n",pixcolname,status);
  fits_write_col(healpix_map_out_fptr,TLONG,1,1L,1L,mapsize,pixel_nb,&status);
  printf("Filled column %s STATUS:%d\n",pixcolname,status);

  fits_insert_col(healpix_map_out_fptr,1,ngalcolname,tform_ngal,&status);
  printf("Inserted column %s STATUS:%d\n",ngalcolname,status);
  fits_write_col(healpix_map_out_fptr,TLONG,1,1L,1L,mapsize,hp_map_value,&status);
  printf("Filled column %s STATUS:%d\n",ngalcolname,status);

  // fits_close_file(catalog_fptr,&status);
  // printf("Close catalog fits file, status %d\n",&status);  

  return 0;
}
