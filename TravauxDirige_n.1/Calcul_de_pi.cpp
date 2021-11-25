#include <mpi.h>
#include <cstdlib>
#include <chrono>
#include <random>

// Attention , ne marche qu'en C++ 11 ou supérieur :
double approximate_pi(unsigned long nbSamples) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution <double > distribution(-1.0, 1.0);
    unsigned long nbDarts = 0;
    // Throw nbSamples darts in the unit square [ -1:1] x [ -1:1]
    for (unsigned sample = 0; sample < nbSamples; ++sample) {
        double x = distribution(generator);
        double y = distribution(generator);
        // Test if the dart is in the unit disk
        if (x * x + y * y <= 1) nbDarts++;
    }
    // Number of nbDarts throwed in the unit disk
    double ratio = double(nbDarts) / double(nbSamples);
    return ratio;
}


int main(int argc, char* argv[]) {
    int nbp;
    int rank;
    int tag = 1234;
    MPI_Status status;

    int total_samples = 99999999;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nbp);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    double r;
    double buf_r;

    int samples = (int)total_samples / nbp + 1;

    if (rank == 0) {

        r = approximate_pi(samples);
        std::cout << "I'm process n° " << rank << " and the calculated Value is: " << 4*r << std::endl;
        for (int k = 1; k <= nbp - 1; k++) {
            MPI_Recv(&buf_r, 1, MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &status);
            r += buf_r;
        }
        double result = 4*r / nbp;
        std::cout << "I'm process n° " << rank << " and the final calculated Value is: " << result << std::endl;
    }
    else {
        r = approximate_pi(samples);
        std::cout << "I'm process n° " << rank << " and the calculated Value is: " << 4*r << std::endl;
        MPI_Send(&r, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;

}
