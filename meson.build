project('pyoculus', 'c', version : '0.4.01',  meson_version: '>= 0.64.0',
)

add_languages('fortran', native: false)

fortran_compiler = meson.get_compiler('fortran')
message('fortran compiler id is ' + fortran_compiler.get_id())


f90flags = []
if fortran_compiler.get_id() == 'gcc'
    f90flags += '-ffree-line-length-none'
    message('adding -ffree-line-length-flag')
endif

f90flags += '-O3'
add_global_arguments(f90flags, language: 'fortran')
    

py3 = import('python').find_installation(pure: false)
py_dep = py3.dependency()

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

sourcefiles = files(
        'pyoculus/problems/SPECfortran/pyvariables.f90',
        'pyoculus/problems/SPECfortran/pybasefn.f90',
        'pyoculus/problems/SPECfortran/pycoords.f90',
        'pyoculus/problems/SPECfortran/pybfield.f90',
        'pyoculus/problems/SPECfortran/pyPJH.f90'
)

module_name = 'pyoculus_spec_fortran_module'


pyoculus_spec_fortran_module_source = custom_target('pyoculus_spec_fortran_module.c',
  input : sourcefiles,  
  output : [module_name+'module.c', module_name+'-f2pywrappers2.f90'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', module_name, '--lower']
)

py3.extension_module(module_name, 
  [sourcefiles, pyoculus_spec_fortran_module_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true
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
