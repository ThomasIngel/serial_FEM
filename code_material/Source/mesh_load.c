#include "hpc.h"

/* load a (index) matrix from a file */
index *mesh_load_index (char *fname, index cols, index *rows)
{
  FILE *file;
  index cnt, j, a, *data;
  
  file = fopen(fname,"r");
  if (file == NULL) return (NULL) ;
  cnt = 0;
  while (fscanf(file,"%zu",&a) != EOF) cnt++;
  rows[0] = cnt/cols;
  if (  cnt - rows[0] * cols ){
    fclose(file);
    printf("\nmesh_load_index() Error!!! cnt = %g, rows %g\n\n", 
            (double) cnt, (double) rows[0]);
    return (NULL);
  }

  data = malloc(cnt * sizeof(index));
  if (!data) return (NULL) ;
  
  fseek(file,0L,SEEK_SET);
  for (j=0; j<cnt; j++) {
    fscanf(file,"%zu",&(data[j]));
  }
  fclose(file);
  return (data);
}

/* load a (double) matrix from a file */
double *mesh_load_double (char *fname, index cols, index *rows)
{
  FILE *file;
  index cnt, j;
  double a, *data;
          
  file = fopen(fname,"r");
  if (file == NULL) return (NULL) ;
  cnt = 0;
  while (fscanf(file,"%lg",&a) != EOF) cnt++;

  rows[0] = cnt/cols;
  if (  cnt - rows[0] * cols ){
    fclose(file);
    printf("\nmesh_load_double() Error!!! cnt = %g, rows %g\n\n", 
            (double) cnt, (double) rows[0]);
    return (NULL);
  }

  data = malloc(cnt * sizeof(double));
  if (!data) return (NULL) ;
  
  fseek(file,0L,SEEK_SET);
  for (j=0; j<cnt; j++) {
    fscanf(file,"%lg",&(data[j]));
  }
  fclose(file);
  return (data);
}

/* load a triplet matrix from a file */
mesh *mesh_load (char *fname)
{
  FILE *file;
  char *tmp;
  index cnt, j, a, *data;
  mesh *M;
  char buffer[512];
  
  M = mesh_alloc(0,0,0) ;               
  if (!M) return (NULL) ;
  M->nedges = 0;

 
  sprintf(buffer,"%s.co",fname);
  printf("Load coordinates from %s\n",buffer);
  M->coord = mesh_load_double(buffer, 2, &(M->ncoord));
  
  sprintf(buffer,"%s.el",fname);
  printf("Load elements from %s\n",buffer);
  M->elem = mesh_load_index(buffer, 7, &(M->nelem));
  
  sprintf(buffer,"%s.bd",fname);
  printf("Load boundary data from %s\n",buffer);
  M->bdry = mesh_load_index(buffer, 4, &(M->nbdry));
  
  return ((!M->coord || !M->elem || !M->bdry) ? mesh_free (M) : M) ;
}
