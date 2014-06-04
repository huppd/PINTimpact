#include "mlguide.hpp"


int main(int argc, char *argv[]){

  int i, N_grids = 20;
  double sol[129], rhs[129];

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
  for (i = 0; i < 129; i++) sol[i] = 0.;
  for (i = 0; i < 129; i++) rhs[i] = 2.;

  //  auto myML = Teuchos::rcp(new bla::MyML<int>(N_grids) );
  auto myML = bla::createMyML<double,int>(N_grids) ;

  myML->apply(sol,rhs);


  /******** End code to set a user-defined smoother ******/
  printf("answer is %e %e %e %e %e\n",sol[0],sol[1],sol[2],sol[3],sol[4]);
  myML = Teuchos::null;


#ifdef ML_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}

