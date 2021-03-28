#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _ChR2H134R_reg(void);
extern void _WB_reg(void);
extern void _WBCN_reg(void);
extern void _WBS_reg(void);
extern void _hhCN_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"mods/ChR2H134R.mod\"");
    fprintf(stderr," \"mods/WB.mod\"");
    fprintf(stderr," \"mods/WBCN.mod\"");
    fprintf(stderr," \"mods/WBS.mod\"");
    fprintf(stderr," \"mods/hhCN.mod\"");
    fprintf(stderr, "\n");
  }
  _ChR2H134R_reg();
  _WB_reg();
  _WBCN_reg();
  _WBS_reg();
  _hhCN_reg();
}
