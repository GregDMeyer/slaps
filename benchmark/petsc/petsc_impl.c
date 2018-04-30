
static char help[] = "Generate a square matrix of given sparsity and size; time MatVec.\n"
                     "Options:\n"
                     "-d\tThe integer dimension of the matrix.\n"
                     "-s\tSparsity (dim/nonzeros) of the matrix.\n"
                     "-i\tThe number of iterations to run MatVec.\n"
                     "-q\tQuiet mode (only output timing).\n";

#include <petscmat.h>
#include <petsctime.h>

int main(int argc,char **args)
{
  Mat             A;                      /* matrix */
  Vec             x, y;
  PetscErrorCode  ierr;
  PetscInt        dim = 100, sparsity = 10, iterations = 100, tot_nz, mpi_size;
  PetscInt        i, j, istart, iend;
  PetscBool       flg, quiet = PETSC_FALSE;
  PetscLogDouble  time;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);
  if (ierr) return ierr;

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  ierr = PetscOptionsGetInt(NULL, NULL, "-d", &dim, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-s", &sparsity, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-i", &iterations, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-q", &quiet, &flg);CHKERRQ(ierr);

  if (!quiet) {
    ierr = PetscPrintf(MPI_COMM_WORLD, "Timing PETSc MatVec.\n");CHKERRQ(ierr);
    ierr = PetscPrintf(MPI_COMM_WORLD, " dim = %d\n sparsity = %d\n iterations = %d\n",
                         dim,       sparsity,       iterations);CHKERRQ(ierr);
  }

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, dim, dim);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);

  /* set 1's in the correct sparsity pattern in the matrix */
  tot_nz = dim/sparsity + 1;
  ierr = MatSeqAIJSetPreallocation(A, tot_nz, NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A, tot_nz/mpi_size + 1, NULL,
                                      tot_nz - tot_nz/mpi_size + 1, NULL);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A, &istart, &iend);CHKERRQ(ierr);

  for (i = istart; i < iend; ++i) {
    /* start from a different spot each time */
    for (j = (91*i) % sparsity; j < dim; j += sparsity) {
      MatSetValue(A, i, j, 1, INSERT_VALUES);
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* now create our vectors */
  ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &y);CHKERRQ(ierr);
  ierr = VecSetSizes(x, PETSC_DECIDE, dim);CHKERRQ(ierr);
  ierr = VecSetSizes(y, PETSC_DECIDE, dim);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecSetFromOptions(y);CHKERRQ(ierr);

  ierr = VecSet(x, 1);CHKERRQ(ierr);

  /* time MatVec */
  ierr = PetscTime(&time);CHKERRQ(ierr);
  for (i = 0; i < iterations; ++i) {
    ierr = MatMult(A, x, y);
  }
  ierr = PetscTimeSubtract(&time);CHKERRQ(ierr);

  if (!quiet) {
    ierr = PetscPrintf(MPI_COMM_WORLD, "Time: ");CHKERRQ(ierr);
  }
  ierr = PetscPrintf(MPI_COMM_WORLD, "%f\n", -time);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}
