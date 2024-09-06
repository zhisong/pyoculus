project('pyoculus', 'c', version : '0.4.01',  meson_version: '>= 0.64.0',
)

add_languages('fortran', native: false)

fortran_compiler = meson.get_compiler('fortran')


f90flags = []
if fortran_compiler.get_id() == 'gcc'
    f90flags += '-ffree-line-length-none'
endif

f90flags += '-O3'
add_project_arguments(f90flags, language: 'fortran')

py3 = import('python').find_installation(pure: false)

sourcefiles = files(
        'pyoculus/problems/SPECfortran/pyvariables.f90',
        'pyoculus/problems/SPECfortran/pybasefn.f90',
        'pyoculus/problems/SPECfortran/pycoords.f90',
        'pyoculus/problems/SPECfortran/pybfield.f90',
        'pyoculus/problems/SPECfortran/pyPJH.f90'
)

python_version = py3.version().split('.')
python_major = python_version[0]
python_minor = python_version[1]
architecture = host_machine.cpu_family()
system = host_machine.system()

# Construct the output filename
### assuming gnu compiler, this gets too unmaintaiable if I make frankenlogic to get the filename
module_name = 'pyoculus_spec_fortran_module'
output_ext = 'cpython-' + python_major + python_minor + '-' + architecture + '-' + system + '-gnu.so'


fortran_sharedobject = custom_target(module_name,
  input : sourcefiles,
  output : module_name +'.'+ output_ext,
  command : [py3, '-m', 'numpy.f2py', '-c', '-m', 'pyoculus_spec_fortran_module', '@INPUT@', '--opt=' + ' '.join(f90flags)],
  install : true,
  install_dir : py3.get_install_dir()
)

py3.install_sources(
           'pyoculus/__init__.py',
            subdir: 'pyoculus' 
            )
py3.install_sources(
           'pyoculus/problems/m3dc1_bfield.py',
           'pyoculus/problems/cylindrical_problem.py',
           'pyoculus/problems/interpolate_coordinates.py',
           'pyoculus/problems/two_waves.py',
           'pyoculus/problems/spec_pjh.py',
           'pyoculus/problems/spec_problem.py',
           'pyoculus/problems/__init__.py',
           'pyoculus/problems/spec_bfield.py',
           'pyoculus/problems/base_problem.py',
           'pyoculus/problems/toroidal_bfield.py',
           'pyoculus/problems/qfm_bfield.py',
           'pyoculus/problems/simsgeo_biot_savart.py',
           'pyoculus/problems/cylindrical_bfield.py',
           'pyoculus/problems/bfield_problem.py',
           'pyoculus/problems/toroidal_problem.py',
           'pyoculus/problems/cartesian_bfield.py',
            subdir: 'pyoculus/problems'
            )
py3.install_sources(
           'pyoculus/integrators/sdc_integrator.py',
           'pyoculus/integrators/base_integrator.py',
           'pyoculus/integrators/__init__.py',
           'pyoculus/integrators/rk_integrator.py',
            subdir: 'pyoculus/integrators' 
            )
py3.install_sources(
           'pyoculus/solvers/__init__.py',
           'pyoculus/solvers/flux_surface_gr.py',
           'pyoculus/solvers/manifolds.py',
           'pyoculus/solvers/base_solver.py',
           'pyoculus/solvers/qfm.py',
           'pyoculus/solvers/lyapunov_exponent.py',
           'pyoculus/solvers/poincare_plot.py',
           'pyoculus/solvers/fixed_point.py',
            subdir: 'pyoculus/solvers' 
            )
py3.install_sources(
           'pyoculus/irrationals/continued_fraction.py',
           'pyoculus/irrationals/__init__.py',
            subdir: 'pyoculus/irrationals' 
)
