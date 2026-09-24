#include "subcell.hpp"
#include "mito.hpp"

void CSubcell::compute_J_ATP_D(void)
{
  // Matches Song's SpatialCell.cu diffusion kernel exactly:
  //   idl/idr/idu/idd/idi/ido computed per-CRU with no-flux boundary
  //   condition idl = (idx==0) ? id : id-1  (missing neighbor -> self),
  //   NOT a mirrored/doubled interior neighbor.
  //   xiATPdiffu = (atpCytp[idl]+atpCytp[idr]-2*atpCytp[id])*D_atp/l_L/l_L
  //              + (atpCytp[idu]+atpCytp[idd]+atpCytp[idi]+atpCytp[ido]-4*atpCytp[id])*D_atp/l_T/l_T

#pragma omp parallel for collapse(3)
  for (int k = 0; k < nz; ++k)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        int id = i + j * nx + k * (nx * ny);

        int idl = (i == 0)      ? id : id - 1;
        int idr = (i == nx - 1) ? id : id + 1;
        int idu = (j == 0)      ? id : id - nx;
        int idd = (j == ny - 1) ? id : id + nx;
        int idi = (k == 0)      ? id : id - nx * ny;
        int ido = (k == nz - 1) ? id : id + nx * ny;

        J_ATP_D[id] = (ATP_cyto[idl] + ATP_cyto[idr] - 2 * ATP_cyto[id]) * DATP / (l_L_atp * l_L_atp)
                    + (ATP_cyto[idu] + ATP_cyto[idd] + ATP_cyto[idi] + ATP_cyto[ido] - 4 * ATP_cyto[id]) * DATP / (l_T_atp * l_T_atp);
      }
    }
  }
}