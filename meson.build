project('pyoculus', 'c', version : '0.4.01',  meson_version: '>= 0.64.0',
)

add_languages('fortran', native: false)

fortran_compiler = meson.get_compiler('fortran')

add_project_arguments('-DFFLAGS=-ffree-line-length-none', language : 'c')

f90flags = []
if fortran_compiler.get_id() == 'gcc'
    f90flags += '-ffree-line-length-none'
endif

f90flags += '-O3'


py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()

incdir_numpy = run_command(
    py3, '-c',
    'import numpy; print(numpy.get_include())',
    check: true,
).stdout().strip()

inc_np = include_directories(incdir_numpy)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)

sourcefiles = files(
        'pyoculus/problems/SPECfortran/pyvariables.f90',
        'pyoculus/problems/SPECfortran/pybasefn.f90',
        'pyoculus/problems/SPECfortran/pycoords.f90',
        'pyoculus/problems/SPECfortran/pybfield.f90',
        'pyoculus/problems/SPECfortran/pyPJH.f90'
)

compiled_fortran = custom_target('pyoculus_spec_fortran_module',
  input : sourcefiles,
  output: 'pyoculus_spec_fortran_modulemodule.c',
  command : [py3, '-m', 'numpy.f2py',  '-m', 'pyoculus_spec_fortran_module', '@INPUT@'],
  install : true,
  install_dir : py3.get_install_dir()
)

py3.extension_module('pyoculus_spec_fortran_module', 
                    sourcefiles, compiled_fortran, 
                    dependencies: fortranobject_dep, 
                    install: true, 
                    fortran_args: f90flags)

py3.install_sources(
           'pyoculus/__init__.py',
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
           'pyoculus/integrators/sdc_integrator.py',
           'pyoculus/integrators/base_integrator.py',
           'pyoculus/integrators/__init__.py',
           'pyoculus/integrators/rk_integrator.py',
           'pyoculus/solvers/__init__.py',
           'pyoculus/solvers/flux_surface_gr.py',
           'pyoculus/solvers/manifolds.py',
           'pyoculus/solvers/base_solver.py',
           'pyoculus/solvers/qfm.py',
           'pyoculus/solvers/lyapunov_exponent.py',
           'pyoculus/solvers/poincare_plot.py',
           'pyoculus/solvers/fixed_point.py',
           'pyoculus/irrationals/continued_fraction.py',
           'pyoculus/irrationals/__init__.py',
            subdir: 'pyoculus' 
)

