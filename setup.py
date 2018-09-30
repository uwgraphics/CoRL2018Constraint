from setuptools import setup
# from geometric_constraints import point_contact_constraint

# pc = point_contact_constraint()
# pc.build_module('point_contact_constraint', 'point_contact_constraint')


setup(name='constraint_recognition',
      version='0.1',
      description='Modeling and recognizing geometric constraints in human demonstration',
      url='',
      author='Guru Subramani',
      author_email='gsubramani@wisc.edu',
      license='MIT',
      packages=['point_contact_constraint','point_on_plane_constraint','prismatic_constraint','planar_constraint',
                'concentric_cylinder_constraint','axial_rotation_constraint'],
      zip_safe=False)