#include "subcell.hpp"
#include "mito.hpp"

void CSubcell::compute_J_ATP_D(void)
{

  // Corners (no-flux boundary: reflected neighbor used in place of the
  // off-grid neighbor, i.e. the same interior neighbor counted twice)
  // x term -> DATP/(l_L_atp^2) (longitudinal), y and z terms -> DATP/(l_T_atp^2) (transverse)
  // written out exactly as Song's formula: D_atp / l_L / l_L  and  D_atp / l_T / l_T
  J_ATP_D[0+0*nx+0*(nx*ny)]=(ATP_cyto[0+0*nx+0*(nx*ny)+1]+ATP_cyto[0+0*nx+0*(nx*ny)+1]-2*ATP_cyto[0+0*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[0+(0+1)*nx+0*(nx*ny)]+ATP_cyto[0+(0+1)*nx+0*(nx*ny)]-2*ATP_cyto[0+0*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[0+0*nx+(0+1)*(nx*ny)]+ATP_cyto[0+0*nx+(0+1)*(nx*ny)]-2*ATP_cyto[0+0*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[0+(ny-1)*nx+(nz-1)*(nx*ny)]=(ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)+1]-2*ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[0+((ny-1)-1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[0+((ny-1)-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[0+(ny-1)*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[0+(ny-1)*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[0+(ny-1)*nx+0*(nx*ny)]=(ATP_cyto[0+(ny-1)*nx+0*(nx*ny)+1]+ATP_cyto[0+(ny-1)*nx+0*(nx*ny)+1]-2*ATP_cyto[0+(ny-1)*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[0+((ny-1)-1)*nx+0*(nx*ny)]+ATP_cyto[0+((ny-1)-1)*nx+0*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[0+(ny-1)*nx+(0+1)*(nx*ny)]+ATP_cyto[0+(ny-1)*nx+(0+1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[0+0*nx+(nz-1)*(nx*ny)]=(ATP_cyto[0+0*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[0+0*nx+(nz-1)*(nx*ny)+1]-2*ATP_cyto[0+0*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[0+(0+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[0+(0+1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[0+0*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[0+0*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[0+0*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[0+0*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[(nx-1)+(ny-1)*nx+0*(nx*ny)]=(ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)-1]+ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)-1]-2*ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[(nx-1)+((ny-1)-1)*nx+0*(nx*ny)]+ATP_cyto[(nx-1)+((ny-1)-1)*nx+0*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[(nx-1)+(ny-1)*nx+(0+1)*(nx*ny)]+ATP_cyto[(nx-1)+(ny-1)*nx+(0+1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[(nx-1)+0*nx+(nz-1)*(nx*ny)]=(ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)-1]+ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[(nx-1)+(0+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[(nx-1)+(0+1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[(nx-1)+0*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[(nx-1)+0*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[(nx-1)+0*nx+0*(nx*ny)]=(ATP_cyto[(nx-1)+0*nx+0*(nx*ny)-1]+ATP_cyto[(nx-1)+0*nx+0*(nx*ny)-1]-2*ATP_cyto[(nx-1)+0*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[(nx-1)+(0+1)*nx+0*(nx*ny)]+ATP_cyto[(nx-1)+(0+1)*nx+0*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[(nx-1)+0*nx+(0+1)*(nx*ny)]+ATP_cyto[(nx-1)+0*nx+(0+1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  J_ATP_D[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)]=(ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)-1]+ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
    (ATP_cyto[(nx-1)+((ny-1)-1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[(nx-1)+((ny-1)-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
    (ATP_cyto[(nx-1)+(ny-1)*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[(nx-1)+(ny-1)*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);

  // x fixed (varies j, then k)
#pragma omp parallel for
  for (int j=1;j<ny-1;j++)
  {
    J_ATP_D[0+j*nx+0*(nx*ny)]=(ATP_cyto[0+j*nx+0*(nx*ny)+1]+ATP_cyto[0+j*nx+0*(nx*ny)+1]-2*ATP_cyto[0+j*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[0+(j+1)*nx+0*(nx*ny)]+ATP_cyto[0+(j-1)*nx+0*(nx*ny)]-2*ATP_cyto[0+j*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[0+j*nx+(0+1)*(nx*ny)]+ATP_cyto[0+j*nx+(0+1)*(nx*ny)]-2*ATP_cyto[0+j*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[(nx-1)+j*nx+0*(nx*ny)]=(ATP_cyto[(nx-1)+j*nx+0*(nx*ny)-1]+ATP_cyto[(nx-1)+j*nx+0*(nx*ny)-1]-2*ATP_cyto[(nx-1)+j*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[(nx-1)+(j+1)*nx+0*(nx*ny)]+ATP_cyto[(nx-1)+(j-1)*nx+0*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[(nx-1)+j*nx+(0+1)*(nx*ny)]+ATP_cyto[(nx-1)+j*nx+(0+1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[0+j*nx+(nz-1)*(nx*ny)]=(ATP_cyto[0+j*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[0+j*nx+(nz-1)*(nx*ny)+1]-2*ATP_cyto[0+j*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[0+(j+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[0+(j-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[0+j*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[0+j*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[0+j*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[0+j*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[(nx-1)+j*nx+(nz-1)*(nx*ny)]=(ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)-1]+ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[(nx-1)+(j+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[(nx-1)+(j-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[(nx-1)+j*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[(nx-1)+j*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
#pragma ivdep
#pragma vector always
    for (int k=1;k<nz-1;k++)
    {
      J_ATP_D[0+j*nx+k*(nx*ny)]=(ATP_cyto[0+j*nx+k*(nx*ny)+1]+ATP_cyto[0+j*nx+k*(nx*ny)+1]-2*ATP_cyto[0+j*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
        (ATP_cyto[0+(j+1)*nx+k*(nx*ny)]+ATP_cyto[0+(j-1)*nx+k*(nx*ny)]-2*ATP_cyto[0+j*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
        (ATP_cyto[0+j*nx+(k+1)*(nx*ny)]+ATP_cyto[0+j*nx+(k-1)*(nx*ny)]-2*ATP_cyto[0+j*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
      J_ATP_D[(nx-1)+j*nx+k*(nx*ny)]=(ATP_cyto[(nx-1)+j*nx+k*(nx*ny)-1]+ATP_cyto[(nx-1)+j*nx+k*(nx*ny)-1]-2*ATP_cyto[(nx-1)+j*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
        (ATP_cyto[(nx-1)+(j+1)*nx+k*(nx*ny)]+ATP_cyto[(nx-1)+(j-1)*nx+k*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
        (ATP_cyto[(nx-1)+j*nx+(k+1)*(nx*ny)]+ATP_cyto[(nx-1)+j*nx+(k-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    }
  }

  // y fixed (varies i, then k, then j)
#pragma omp parallel for
  for (int i=1;i<(nx-1);i++)
  {
    J_ATP_D[i+0*nx+0*(nx*ny)]=(ATP_cyto[i+0*nx+0*(nx*ny)+1]+ATP_cyto[i+0*nx+0*(nx*ny)-1]-2*ATP_cyto[i+0*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[i+(0+1)*nx+0*(nx*ny)]+ATP_cyto[i+(0+1)*nx+0*(nx*ny)]-2*ATP_cyto[i+0*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[i+0*nx+(0+1)*(nx*ny)]+ATP_cyto[i+0*nx+(0+1)*(nx*ny)]-2*ATP_cyto[i+0*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[i+(ny-1)*nx+0*(nx*ny)]=(ATP_cyto[i+(ny-1)*nx+0*(nx*ny)+1]+ATP_cyto[i+(ny-1)*nx+0*(nx*ny)-1]-2*ATP_cyto[i+(ny-1)*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[i+((ny-1)-1)*nx+0*(nx*ny)]+ATP_cyto[i+((ny-1)-1)*nx+0*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[i+(ny-1)*nx+(0+1)*(nx*ny)]+ATP_cyto[i+(ny-1)*nx+(0+1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[i+0*nx+(nz-1)*(nx*ny)]=(ATP_cyto[i+0*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[i+0*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[i+0*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[i+(0+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[i+(0+1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[i+0*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[i+0*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[i+0*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[i+0*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[i+(ny-1)*nx+(nz-1)*(nx*ny)]=(ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[i+((ny-1)-1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[i+((ny-1)-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[i+(ny-1)*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[i+(ny-1)*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
#pragma ivdep
#pragma vector always
    for (int k=1;k<nz-1;k++)
    {
      J_ATP_D[i+0*nx+k*(nx*ny)]=(ATP_cyto[i+0*nx+k*(nx*ny)+1]+ATP_cyto[i+0*nx+k*(nx*ny)-1]-2*ATP_cyto[i+0*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
        (ATP_cyto[i+(0+1)*nx+k*(nx*ny)]+ATP_cyto[i+(0+1)*nx+k*(nx*ny)]-2*ATP_cyto[i+0*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
        (ATP_cyto[i+0*nx+(k+1)*(nx*ny)]+ATP_cyto[i+0*nx+(k-1)*(nx*ny)]-2*ATP_cyto[i+0*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
      J_ATP_D[i+(ny-1)*nx+k*(nx*ny)]=(ATP_cyto[i+(ny-1)*nx+k*(nx*ny)+1]+ATP_cyto[i+(ny-1)*nx+k*(nx*ny)-1]-2*ATP_cyto[i+(ny-1)*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
        (ATP_cyto[i+((ny-1)-1)*nx+k*(nx*ny)]+ATP_cyto[i+((ny-1)-1)*nx+k*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
        (ATP_cyto[i+(ny-1)*nx+(k+1)*(nx*ny)]+ATP_cyto[i+(ny-1)*nx+(k-1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    }
    // z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<ny-1;j++)
    {
      J_ATP_D[i+j*nx+0*(nx*ny)]=(ATP_cyto[i+j*nx+0*(nx*ny)+1]+ATP_cyto[i+j*nx+0*(nx*ny)-1]-2*ATP_cyto[i+j*nx+0*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
        (ATP_cyto[i+(j+1)*nx+0*(nx*ny)]+ATP_cyto[i+(j-1)*nx+0*(nx*ny)]-2*ATP_cyto[i+j*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
        (ATP_cyto[i+j*nx+(0+1)*(nx*ny)]+ATP_cyto[i+j*nx+(0+1)*(nx*ny)]-2*ATP_cyto[i+j*nx+0*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
      J_ATP_D[i+j*nx+(nz-1)*(nx*ny)]=(ATP_cyto[i+j*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[i+j*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[i+j*nx+(nz-1)*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
        (ATP_cyto[i+(j+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[i+(j-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[i+j*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
        (ATP_cyto[i+j*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[i+j*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[i+j*nx+(nz-1)*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    }
  }

  // remaining k-fixed edges
#pragma ivdep
#pragma vector always
  for (int k=1;k<nz-1;k++)
  {
    J_ATP_D[0+0*nx+k*(nx*ny)]=(ATP_cyto[0+0*nx+k*(nx*ny)+1]+ATP_cyto[0+0*nx+k*(nx*ny)+1]-2*ATP_cyto[0+0*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[0+(0+1)*nx+k*(nx*ny)]+ATP_cyto[0+(0+1)*nx+k*(nx*ny)]-2*ATP_cyto[0+0*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[0+0*nx+(k+1)*(nx*ny)]+ATP_cyto[0+0*nx+(k-1)*(nx*ny)]-2*ATP_cyto[0+0*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[0+(ny-1)*nx+k*(nx*ny)]=(ATP_cyto[0+(ny-1)*nx+k*(nx*ny)+1]+ATP_cyto[0+(ny-1)*nx+k*(nx*ny)+1]-2*ATP_cyto[0+(ny-1)*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[0+((ny-1)-1)*nx+k*(nx*ny)]+ATP_cyto[0+((ny-1)-1)*nx+k*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[0+(ny-1)*nx+(k+1)*(nx*ny)]+ATP_cyto[0+(ny-1)*nx+(k-1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[(nx-1)+0*nx+k*(nx*ny)]=(ATP_cyto[(nx-1)+0*nx+k*(nx*ny)-1]+ATP_cyto[(nx-1)+0*nx+k*(nx*ny)-1]-2*ATP_cyto[(nx-1)+0*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[(nx-1)+(0+1)*nx+k*(nx*ny)]+ATP_cyto[(nx-1)+(0+1)*nx+k*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[(nx-1)+0*nx+(k+1)*(nx*ny)]+ATP_cyto[(nx-1)+0*nx+(k-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
    J_ATP_D[(nx-1)+(ny-1)*nx+k*(nx*ny)]=(ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)-1]+ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)-1]-2*ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
      (ATP_cyto[(nx-1)+((ny-1)-1)*nx+k*(nx*ny)]+ATP_cyto[(nx-1)+((ny-1)-1)*nx+k*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
      (ATP_cyto[(nx-1)+(ny-1)*nx+(k+1)*(nx*ny)]+ATP_cyto[(nx-1)+(ny-1)*nx+(k-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
  }

  // interior — covers the vast majority of the grid
#pragma omp parallel for
  for (int k=1;k<nz-1;k++)
  {
    for (int j=1;j<ny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nx-1;i++)
      {
        J_ATP_D[i+j*nx+k*(nx*ny)]=(ATP_cyto[i+j*nx+k*(nx*ny)+1]+ATP_cyto[i+j*nx+k*(nx*ny)-1]-2*ATP_cyto[i+j*nx+k*(nx*ny)])*DATP/(l_L_atp*l_L_atp)+
          (ATP_cyto[i+(j+1)*nx+k*(nx*ny)]+ATP_cyto[i+(j-1)*nx+k*(nx*ny)]-2*ATP_cyto[i+j*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp)+
          (ATP_cyto[i+j*nx+(k+1)*(nx*ny)]+ATP_cyto[i+j*nx+(k-1)*(nx*ny)]-2*ATP_cyto[i+j*nx+k*(nx*ny)])*DATP/(l_T_atp*l_T_atp);
      }
    }
  }
}