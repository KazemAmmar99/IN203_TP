#include <mpi.h>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
	time_t t;

	srand((unsigned) time(&t));
	int nbp;
	int rank;
	int tag = 1234;

	MPI_Status status;

	int token = rand() % 100 + 1;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &nbp);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		std::cout << "I'm process n°" << rank << " and the message value is : " << token << std::endl;
		token += 1;
		MPI_Send(&token, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
		MPI_Recv(&token, 1, MPI_INT, nbp-1, tag, MPI_COMM_WORLD, &status);
		std::cout << "I'm process n°" << rank << " and the final value is : " << token << std::endl;

	}
	else {
		MPI_Recv(&token, 1, MPI_INT, rank - 1, tag, MPI_COMM_WORLD, &status);
		std::cout << "I'm process n°" << rank << " and the message value is : " << token << std::endl;
		token += 1;
		MPI_Send(&token, 1, MPI_INT, (rank+1)%nbp, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;

}