#include "docopt.h"
#include <iostream>

static const char USAGE[] =
R"(

PURIFY is a software package that applies state-of-the-art sparse reconstruction algorithms
to radio interferometric measurements. The algorithms use wavelet dictionaries to construct
an estimate of the sky brightness distribution, and can be distributed over clusters for
large data sets.
Usage:
  purify <input_name> <output_name> [options]

  <input_name>  Name of input measurement set or visibility file.
  <output_name>  Prefix of output images.

Options:
  -h --help     Show this screen.
  --version     Show version.
  --niters <iters>     Number of iterations [default: 100].
  --logging_level <logging>  Level of logging output [default: debug].
  --weighting <weighting>     Type of weighting scheme to use [default: natural].
  --stokes <stokes>     Choice of stokes parameter to image [default: I].
  --noisefile <noisefile>     File to estimate noise using Stokes V [default: ""].
  --channel_averaging <channel_averaging>     The number of channels to average over [default: 0].
  --beta <beta>     The scaling of the step size gamma [default: 1e-3].
  --over_sample <over_sample_factor>     The over sampling factor for the FFT in gridding/degridding [default: 2].
  --kernel <kernel>     The choice of interpolation kernel [default: kb].
  --kernel_support <J>     The size of the interpolation kernel support [default: 4].
  --width <width>     Image width in pixels [default: 512].
  --height <height>     Image height in pixels [default: 512].
  --cellsizex <cellsizex>     Pixel size in arcseconds [default: 0].
  --cellsizey <cellsizey>     Pixel size in arcseconds [default: 0].
  --primary_beam <primary_beam>     Add primary beam of a telescope in the measurement operator [default: none].
  --fft_grid_correction     Use FFT to calculate the interpolation kernel correction [default: false].
  --fftw_plan <plan>     Type of planning for FFTW [default: measure].
  --gradient <dim>     Add gradient operator along x or y axis to measurement operator [default: none].
  --use_w_term     Use w-projection method [default: false].
  --energy_fraction_chirp <energy_fraction_chirp>     Fraction of clipping of the chirp in w-projection method [default: 0.99999].
  --energy_fraction_wproj <energy_fraction_wproj>     Fraction of clipping of the energy in interpolation kernel with w-projection [default: 0.99999].
  --algo_update     Use lambda function to record or update algorithm variables [default: true].
  --update_output     Save output after each iteration [default: false].
  --adapt_gamma     Update gamma/stepsize [default: true].
  --run_diagnostic  Save and output diagnostic information [default: false].
  --warmstart     If to use warmstart, not sure if this would work with reweighting as implimented [default: false].
  --no_reweighted     If to use reweighting [default: true].
  --relative_gamma_adapt <rel_diff>     When the relative difference of a step size update is more than this value, the step size will accept the update [default: 0.01].
  --adapt_iter <adapt_iters>     The number of iterations to adapt the step size. The step size will remain constant above this limit of iterations [default: 100].
  --power_method_iterations <power_iter>     The limit of iterations in the power method used to normalise the measurement operator [default: 100].
  --l2_bound <l2_bound>     Factor to multiply scale the l2 bound (epsilon) by [default: 1.4].
  --relative_variation <relative_variation>     Relative difference in model for convergence [default: 5e-3].
  --residual_convergence     Factor of epsilon the l2 norm residuals can have for convergence, -1 means it will choose epsilon by default [default: 1].
  --positive     Use positivity constraint [default: true]. # If default is true? how can this be set as false?

)";

int main(int argc, const char** argv)
{
    std::map<std::string, docopt::value> args
        = docopt::docopt(USAGE,
                         { argv + 1, argv + argc },
                         true,               // show help if requested
                         "2.0");  // version string

    for(auto const& arg : args) {
      std::cout << arg.first << " " <<   arg.second << std::endl;
    }

    return 0;
}
