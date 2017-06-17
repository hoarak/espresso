int ObservableComPosition::actual_calculate() {
  double* A = last_value;
  double p_com[3] = { 0. , 0., 0. } ;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    p_com[0] += (partCfg[ids->e[i]]).p.mass*partCfg[ids->e[i]].pos()[0];
    p_com[1] += (partCfg[ids->e[i]]).p.mass*partCfg[ids->e[i]].pos()[1];
    p_com[2] += (partCfg[ids->e[i]]).p.mass*partCfg[ids->e[i]].pos()[2];
    total_mass += (partCfg[ids->e[i]]).p.mass;
  }
  A[0]=p_com[0]/total_mass;
  A[1]=p_com[1]/total_mass;
  A[2]=p_com[2]/total_mass;
  return 0;
}
