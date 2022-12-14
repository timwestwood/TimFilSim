// seeding.cu

#include "seeding.hpp"

#if SURFACE_OF_REVOLUTION_BODIES

  #include <iostream>
  #include <fstream>
  #include <cmath>
  #include <random>
  #include <string>
  #include "matrix.hpp"

  __global__ void find_nearest_neighbours(int *const ids, const double *const samples, const int num_samples, const double *const X, const int N){

    // Work out which sample(s) this thread will compute the nearest-neighbour for
    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    // Declare the shared memory for this thread block
    __shared__ double x_shared[THREADS_PER_BLOCK];
    __shared__ double y_shared[THREADS_PER_BLOCK];
    __shared__ double z_shared[THREADS_PER_BLOCK];

    double xi, yi, zi;
    int min_id;
    double min_dist = 1e100;

    for (int i = index; (i-threadIdx.x) < num_samples; i+=stride){

      if (i < num_samples){

        xi = samples[3*i];
        yi = samples[3*i + 1];
        zi = samples[3*i + 2];

      }

      for (int j_start = 0; j_start < N; j_start += THREADS_PER_BLOCK){

        const int j_to_read = j_start + threadIdx.x;

        if (j_to_read < N){

          x_shared[threadIdx.x] = X[3*j_to_read];
          y_shared[threadIdx.x] = X[3*j_to_read + 1];
          z_shared[threadIdx.x] = X[3*j_to_read + 2];

        }

        __syncthreads();

        if (i < num_samples){

          for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < N); j++){

            const double dx = xi - x_shared[j];
            const double dy = yi - y_shared[j];
            const double dz = zi - z_shared[j];

            const double dist = sqrt(dx*dx + dy*dy + dz*dz);

            if (dist < min_dist){

              min_dist = dist;
              min_id = j_start + j;

            }

          }

        }

        __syncthreads();

      }

      if (i < num_samples){

        ids[i] = min_id;

      }

    }

  } // End of kernel.






  class shape_fourier_description{

  public:

  int num_modes;
  matrix cos_coeffs;
  matrix sin_coeffs;

  std::mt19937 gen;

  matrix polar_angle_area_fractions;

  ~shape_fourier_description(){};

  shape_fourier_description(){

    std::ifstream input_file(GENERATRIX_FILE_NAME+std::string(".fourier_modes"));

    if (input_file.good()){

      input_file >> num_modes;

      cos_coeffs = matrix(1, num_modes);
      sin_coeffs = matrix(1, num_modes);

      for (int n = 0; n < num_modes; n++){

        input_file >> cos_coeffs(n);
        input_file >> sin_coeffs(n);

      }

    } else {

      std::cout << std::endl << std::endl << " The required generatrix file '" << GENERATRIX_FILE_NAME << ".fourier_modes" << "' was not found." << std::endl;
      exit(-1);

    }

    input_file.close();

    if (std::string(GENERATRIX_FILE_NAME) != std::string("sphere")){

      input_file.open(GENERATRIX_FILE_NAME + std::string(".polar_angle_area_fractions"));

      if (input_file.good()){

        int num_theta_locs;

        input_file >> num_theta_locs;

        polar_angle_area_fractions = matrix(1, num_theta_locs);

        for (int n = 0; n < num_theta_locs; n++){

          input_file >> polar_angle_area_fractions(n);

        }

      } else {

        find_polar_angle_area_fractions();

      }

      input_file.close();

    }

    // Random seed for the RNG.
    std::random_device rd{};
    gen = std::mt19937{rd()};

  };

  double radius(const double theta) const {

    matrix cos_mat(num_modes,1);
    matrix sin_mat(num_modes,1);

    for (int n = 0; n < num_modes; n++){

      cos_mat(n) = std::cos(n*theta);
      sin_mat(n) = std::sin(n*theta);

    }

    return cos_coeffs*cos_mat + sin_coeffs*sin_mat;

  };

  matrix location(const double theta, const double phi) const {

    matrix loc(3,1);

    const double r = radius(theta);

    const double cp = std::cos(phi);
    const double sp = std::sin(phi);

    const double ct = std::cos(theta);
    const double st = std::sin(theta);

    loc(0) = r*cp*st;
    loc(1) = r*sp*st;
    loc(2) = r*ct;

    return loc;

  };

  matrix azi_dir(const double phi) const {

    matrix e(3,1);

    e(0) = -std::sin(phi);
    e(1) = std::cos(phi);
    e(2) = 0.0;

    return e;

  };

  double radius_deriv(const double theta) const {

    double out = 0.0;

    for (int n = 1; n < num_modes; n++){

        out += n*(sin_coeffs(n)*std::cos(n*theta) - cos_coeffs(n)*std::sin(n*theta));

    }

    return out;

  };

  double radius_second_deriv(const double theta) const {

    double out = 0.0;

    for (int n = 1; n < num_modes; n++){

        out -= n*n*(sin_coeffs(n)*std::sin(n*theta) + cos_coeffs(n)*std::cos(n*theta));

    }

    return out;

  };

  matrix polar_dir(const double theta, const double phi) const {

    const double r = radius(theta);
    const double rprime = radius_deriv(theta);

    const double cp = std::cos(phi);
    const double sp = std::sin(phi);

    const double ct = std::cos(theta);
    const double st = std::sin(theta);

    matrix e(3,1);

    const double temp = rprime*st + r*ct;
    e(0) = cp*temp;
    e(1) = sp*temp;
    e(2) = rprime*ct - r*st;

    return e/norm(e);

  };

  matrix surface_normal(const double theta, const double phi) const {

    return cross(polar_dir(theta, phi), azi_dir(phi));

  };

  matrix full_frame(const double theta, const double phi) const {

    matrix frame(9,1);

    frame.set_block(0, 3, polar_dir(theta, phi));
    frame.set_block(3, 3, azi_dir(phi));
    frame.set_block(6, 3, cross(frame.get_block(0,3), frame.get_block(3,3)));

    return frame;

  };

  matrix random_point(){

    // N.B. The Mersenne twister "gen" is modified whenever we call upon it, so this method cannot be
    // const and consequently this class cannot be passed as const to the seeding function.

    matrix sample(3,1);

    if (std::string(GENERATRIX_FILE_NAME) == std::string("sphere")){

      // There's a really simple way of generating uniform samples on a sphere.
      std::normal_distribution<double> d(0.0, 1.0);

      sample(0) = d(gen);
      sample(1) = d(gen);
      sample(2) = d(gen);

      const double sample_norm = norm(sample);

      sample /= sample_norm;

    } else {

      std::uniform_real_distribution<double> d(0.0, 1.0);

      const double u1 = d(gen);
      const double phi = 2.0*PI*u1;

      const double u2 = d(gen);

      const int num_theta = polar_angle_area_fractions.num_cols * polar_angle_area_fractions.num_rows; // So it doesn't matter whether I made it a row vector or a column vector.

      for (int n = 1; n < num_theta; n++){

        if (u2 < polar_angle_area_fractions(n)){

          const double frac = (u2 - polar_angle_area_fractions(n-1))/(polar_angle_area_fractions(n) - polar_angle_area_fractions(n-1));
          const double theta_l = PI*(n-1)/double(num_theta-1);
          const double dtheta = PI/double(num_theta-1);

          const double theta = theta_l + frac*dtheta;

          sample = location(theta, phi);

          break;

        }

      }

    }

    return sample;

  };

  double polar_angle_area_fraction_integrand(const double theta) const {

    const double r = radius(theta);
    const double rprime = radius_deriv(theta);

    return r*std::sin(theta)*std::sqrt(r*r + rprime*rprime);

  };

  void find_polar_angle_area_fractions(){

    int num_theta_locs = 1000;

    polar_angle_area_fractions = matrix(1, num_theta_locs);

    polar_angle_area_fractions(0) = 0.0;

    double F_old = 0.0;

    for (int n = 1; n < num_theta_locs; n++){

      const double theta = PI*n/double(num_theta_locs-1);
      const double F = polar_angle_area_fraction_integrand(theta);
      polar_angle_area_fractions(n) = polar_angle_area_fractions(n-1) + 0.5*PI*(F + F_old)/double(num_theta_locs-1); // Trapezium rule.

      F_old = F;

    }

    // Scale by the total (which is actually the total area of the shape divided by 2pi)
    for (int n = 1; n < num_theta_locs-1; n++){

      polar_angle_area_fractions(n) /= polar_angle_area_fractions(num_theta_locs-1);

    }

    polar_angle_area_fractions(num_theta_locs-1) = 1.0; // Do this to avoid any rounding-error fishyness.

    // Write to file
    std::ofstream file(GENERATRIX_FILE_NAME+std::string(".polar_angle_area_fractions"));

    file << num_theta_locs << " ";

    for (int n = 0; n < num_theta_locs; n++){

      file << polar_angle_area_fractions(n) << " ";

    }

    file.close();

  };

  double surface_projection_error(const double theta, const double x1, const double x2) const {

    const double r = radius(theta);
    const double rp = radius_deriv(theta);

    const double c = std::cos(theta);
    const double s = std::sin(theta);

    return r*rp + x1*(r*s - rp*c) - x2*(r*c + rp*s);

  };

  double surface_projection_error_deriv(const double theta, const double x1, const double x2) const {

    const double r = radius(theta);
    const double rp = radius_deriv(theta);
    const double rpp = radius_second_deriv(theta);

    const double c = std::cos(theta);
    const double s = std::sin(theta);

    return r*rpp + rp*rp + x1*((r-rpp)*c + 2.0*rp*s) - x2*((rpp-r)*s + 2.0*rp*c);

  };

  void project_onto_surface(double *const X) const {

    double theta = std::atan2(std::sqrt(X[0]*X[0] + X[1]*X[1]), X[2]);

    if (std::string(GENERATRIX_FILE_NAME) == std::string("sphere")){

      const double curr_radius = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
      const double targ_radius = radius(theta);

      X[0] *= targ_radius/curr_radius;
      X[1] *= targ_radius/curr_radius;
      X[2] *= targ_radius/curr_radius;

    } else {

      const double x1 = X[2];
      const double x2 = std::sqrt(X[0]*X[0] + X[1]*X[1]);

      double error = surface_projection_error(theta, x1, x2);

      while (error > 1e-2){

        const double error_prime = surface_projection_error_deriv(theta, x1, x2);

        theta -= error/error_prime; // Newton's method.

        error = surface_projection_error(theta, x1, x2);

      }

      const double phi = std::atan2(X[1], X[0]);

      const matrix loc = location(theta, phi);

      X[0] = loc(0);
      X[1] = loc(1);
      X[2] = loc(2);

    }

  };

  };







  void equal_area_seeding(double *const pos_ref, double *const polar_dir_refs, double *const azi_dir_refs, double *const normal_refs, const int N, shape_fourier_description& shape){

    if (N == 0){

      return; // Because of the squirmer-style simulations, the code may try to seed 0 filaments.

    }

    // This function essentially implements a version of MacQueen's algorithm.
    // We seed a large number of points randomly on the surface according to a uniform distribution, and then find
    // our N points as the centres of regions of (approximately) equal area through N-means clustering.

    int samples_per_iter = 1000;

    // Allocate memory for the CUDA nearest-neighbour search.
    double *samples, *X;
    cudaMallocManaged(&samples, 3*samples_per_iter*sizeof(double));
    cudaMallocManaged(&X, 3*N*sizeof(double));

    int *sample_nn_ids;
    cudaMallocManaged(&sample_nn_ids, samples_per_iter*sizeof(int));

    const int num_thread_blocks = (samples_per_iter + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    // Generate initial guesses.
    if (std::string(GENERATRIX_FILE_NAME) == std::string("sphere")){

      // We use the spiral distribution given by Saff and Kuijlaars (1997) as our initial positions.
      double phi = 0.0;

      for (int n = 0; n < N; n++){

        const double theta = std::acos((N == 1) ? -1.0 : 2.0*n/double(N-1) - 1.0);

        const matrix d = shape.location(theta, phi);

        X[3*n] = d(0);
        X[3*n + 1] = d(1);
        X[3*n + 2] = d(2);

        if (n == N-2){

          phi = 0.0;

        } else {

          const double r = std::sqrt(d(0)*d(0) + d(1)*d(1));

          phi += 3.6/(r*std::sqrt(N));

        }

      }

    } else {

      // Use random initial positions.
      for (int n = 0; n < N; n++){

        const matrix sample = shape.random_point();
        X[3*n] = sample(0);
        X[3*n + 1] = sample(1);
        X[3*n + 2] = sample(2);

      }

    }

    std::vector<int> count(N);
    std::vector<int> old_count(N); // Used to check for convergence.

    for (int n = 0; n < N; n++){

      count[n] = 1;
      old_count[n] = 1;

    }

    int num_iters = 0;

    while (num_iters < 10000){

      // Sample from the uniform distribution on the surface.
      for (int n = 0; n < samples_per_iter; n++){

        const matrix sample = shape.random_point();

        samples[3*n] = sample(0);
        samples[3*n + 1] = sample(1);
        samples[3*n + 2] = sample(2);

      }
      

      // Find the candidate points closest to the random samples.
      find_nearest_neighbours<<<num_thread_blocks, THREADS_PER_BLOCK>>>(sample_nn_ids, samples, samples_per_iter, X, N);
      cudaDeviceSynchronize();

      // Update the candidate points.
      for (int n = 0; n < samples_per_iter; n++){

        const int min_id = sample_nn_ids[n];

        X[3*min_id] = (count[min_id]*X[3*min_id] + samples[3*n])/double(count[min_id] + 1);
        X[3*min_id + 1] = (count[min_id]*X[3*min_id + 1] + samples[3*n + 1])/double(count[min_id] + 1);
        X[3*min_id + 2] = (count[min_id]*X[3*min_id + 2] + samples[3*n + 2])/double(count[min_id] + 1);

        count[min_id]++;

        shape.project_onto_surface(&X[3*min_id]);

      }

      num_iters++;

      // Check for convergence.
      // In general, I only expect this to trigger an early exit for small numbers of points (e.g. some filament seeding problems),
      // and even then possibly only on spheres (or any other shapes that actually have some exact solutions).
      if (num_iters % 100 == 0){

        int max_change = 0;
        int min_change = 1000000;

        for (int n = 0; n < N; n++){

          if (max_change < count[n]-old_count[n]){

            max_change = count[n] - old_count[n];

          }

          if (min_change > count[n]-old_count[n]){

            min_change = count[n] - old_count[n];

          }

        }

        old_count = count;

        if (double(min_change)/double(max_change) > 0.99){

          break;

        }

      }

    }

    // Write the data for the final positions
    for (int n = 0; n < N; n++){

      pos_ref[3*n] = X[3*n];
      pos_ref[3*n + 1] = X[3*n + 1];
      pos_ref[3*n + 2] = X[3*n + 2];

      const double theta = std::atan2(std::sqrt(X[3*n]*X[3*n] + X[3*n + 1]*X[3*n + 1]), X[3*n + 2]);
      const double phi = std::atan2(X[3*n + 1], X[3*n]);

      matrix frame = shape.full_frame(theta, phi);

      polar_dir_refs[3*n] = frame(0);
      polar_dir_refs[3*n + 1] = frame(1);
      polar_dir_refs[3*n + 2] = frame(2);

      azi_dir_refs[3*n] = frame(3);
      azi_dir_refs[3*n + 1] = frame(4);
      azi_dir_refs[3*n + 2] = frame(5);

      normal_refs[3*n] = frame(6);
      normal_refs[3*n + 1] = frame(7);
      normal_refs[3*n + 2] = frame(8);

    }

  };





  void seed_blobs(double *const blob_references, double *const polar_dir_refs, double *const azi_dir_refs, double *const normal_refs){

    const std::string file_name_trunk = GENERATRIX_FILE_NAME+std::to_string(NBLOB);

    std::cout << std::endl << std::endl << "Input file " << "'" << file_name_trunk << ".seed'" << " was not found or could not be opened." << std::endl;
    std::cout << "Seeking an equal-area distribution for the blobs..." << std::endl;

    shape_fourier_description shape;

    equal_area_seeding(blob_references, polar_dir_refs, azi_dir_refs, normal_refs, NBLOB, shape);

    std::ofstream blob_ref_file(file_name_trunk + ".seed");
    std::ofstream polar_file(file_name_trunk + ".polar_dir");
    std::ofstream azi_file(file_name_trunk + ".azi_dir");
    std::ofstream normal_file(file_name_trunk + ".normal");

    for (int n = 0; n < 3*NBLOB; n++){

      blob_ref_file << blob_references[n] << " ";
      polar_file << polar_dir_refs[n] << " ";
      azi_file << azi_dir_refs[n] << " ";
      normal_file << normal_refs[n] << " ";

    }

    blob_ref_file.close();
    polar_file.close();
    azi_file.close();
    normal_file.close();

    std::cout << "...done!" << std::endl;

  };





  void seed_filaments(double *const filament_references, double *const polar_dir_refs, double *const azi_dir_refs, double *const normal_refs){

    std::string file_name_trunk = GENERATRIX_FILE_NAME+std::to_string(NFIL);

    #if EQUATORIAL_SEEDING

      file_name_trunk += "_equatorial";

    #elif PLATY_SEEDING

      file_name_trunk += "_platy";

    #endif

    std::cout << std::endl << std::endl << "Input file " << "'" << file_name_trunk << ".seed'" << " was not found or could not be opened." << std::endl;

    shape_fourier_description shape;

    #if EQUATORIAL_SEEDING

    // TODO: Make this work for non-spheres.

      std::cout << "Seeding locations according to a spiral pattern..." << std::endl;

      const double theta_max = 0.5/double(R_OVER_L); // = 0.5 * 2*PI*L/(2*PI*R), meaning the width of the band is L as measured across the sphere surface.
      const double h = 2.0*std::sin(theta_max);
      const int num_fils_for_sphere = std::ceil(NFIL*2.0/h);
      const int offset = std::round(0.5*(num_fils_for_sphere - NFIL));

      double phi = 0.0;

      for (int n = 0; n < NFIL; n++){

        filament_references[3*n + 2] = 2.0*(n + offset)/double(num_fils_for_sphere - 1) - 1.0;
        const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);
        filament_references[3*n] = r*std::cos(phi);
        filament_references[3*n + 1] = r*std::sin(phi);

        phi += 3.6/(r*std::sqrt(num_fils_for_sphere));

      }

    #elif PLATY_SEEDING

    // TODO: Make this work for non-spheres.

      std::cout << "Seeding..." << std::endl;

      const double theta_max = 0.5/double(R_OVER_L); // = 0.5 * 2*PI*L/(2*PI*R), meaning the width of the band is L as measured across the sphere surface.
      const double h = 2.0*std::sin(theta_max);
      const int num_fils_in_main_band = std::round(0.9*NFIL); // Put 90% in the main band
      const int num_fils_for_sphere = std::ceil(num_fils_in_main_band*2.0/h);
      const int offset = std::round(0.5*(num_fils_for_sphere - num_fils_in_main_band));

      // Seed the main band
      double phi = 0.0;

      for (int n = 0; n < num_fils_in_main_band; n++){

        filament_references[3*n + 2] = 2.0*(n + offset)/double(num_fils_for_sphere - 1) - 1.0;
        const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);
        filament_references[3*n] = r*std::cos(phi);
        filament_references[3*n + 1] = r*std::sin(phi);

        phi += 3.6/(r*std::sqrt(num_fils_for_sphere));

      }

      // Seed the remaining 10% at the bottom of the swimmer
      for (int n = num_fils_in_main_band; n < NFIL; n++){

        filament_references[3*n + 2] = -0.9; // Some fixed location for the bottom band
        const double r = std::sqrt(1.0 - filament_references[3*n + 2]*filament_references[3*n + 2]);
        phi = 2.0*PI*double(n - num_fils_in_main_band)/double(NFIL - num_fils_in_main_band);
        filament_references[3*n] = r*std::cos(phi);
        filament_references[3*n + 1] = r*std::sin(phi);

      }

    #elif UNIFORM_SEEDING

      std::cout << "Seeking an equal-area distribution for the filaments..." << std::endl;

      equal_area_seeding(filament_references, polar_dir_refs, azi_dir_refs, normal_refs, NFIL, shape);

    #endif

    std::ofstream fil_ref_file(file_name_trunk + ".seed");
    std::ofstream polar_file(file_name_trunk + ".polar_dir");
    std::ofstream azi_file(file_name_trunk + ".azi_dir");
    std::ofstream normal_file(file_name_trunk + ".normal");

    for (int n = 0; n < 3*NFIL; n++){

      fil_ref_file << filament_references[n] << " ";
      polar_file << polar_dir_refs[n] << " ";
      azi_file << azi_dir_refs[n] << " ";
      normal_file << normal_refs[n] << " ";

    }

    fil_ref_file.close();
    polar_file.close();
    azi_file.close();
    normal_file.close();

    std::cout << "...done!" << std::endl;

  };

#endif
