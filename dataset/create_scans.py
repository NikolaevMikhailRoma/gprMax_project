import os
import numpy as np
import matplotlib.pyplot as plt
from gprMax.gprMax import api
from tools.plot_Bscan import mpl_plot
from tools.plot_Ascan import mpl_plot as mpl_Ascan_plot

from tools.outputfiles_merge import get_output_data, merge_files


def material_peplinski_box(f1=0.5, f2=0.5, f3=2.0, f4=2.66,
                           f5=0.001, f6=0.25,
                           material_name=None,
                           random=False):
    """
    # soil_peplinski: f1 f2 f3 f4 f5 f6 str1
    • f1 is the sand fraction of the soil.
    • f2 is the clay fraction of the soil.
    • f3 is the bulk density of the soil in grams per centimetre cubed.
    • f4 is the density of the sand particles in the soil in grams per centimetre cubed.
    • f5 and f6 define a range for the volumetric water fraction of the soil.
    • str1 is an identifier for the soil.
    For example for a soil with sand fraction 0.5, clay fraction 0.5, bulk density 2 𝑔/𝑐𝑚3, sand particle density of
    2.66 𝑔/𝑐𝑚3, and a volumetric water fraction range of 0.001 - 0.25 use:
    #soil_peplinski: 0.5 0.5 2.0 2.66 0.001 0.25 my_soil.
    """
    if material_name == None:
        material_name = '1'
    txt = '#soil_peplinski: {} {} {} {} {} {} {}\n\n'.format(f1, f2, f3, f4, f5, f6, material_name)
    return txt


def geometry_peplinski_box(x1, y1, z1, str1, str2, x0=0, y0=0, z0=0,
                           f7=1.5, f8=1, f9=1, f10=1, i1=50, ):
    """
    #fractal_box: f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 i1 str1 str2 [i2] [c1]

    f1 f2 f3 are the lower left (x,y,z) coordinates of the parallelepiped,
    and f4 f5 f6 are the upper right(x,y,z) coordinates of the parallelepiped.
    f7 is the fractal dimension which, for an orthogonal parallelepiped, should take values between zero and three.
    f8 is used to weight the fractal in the x direction.
    f9 is used to weight the fractal in the y direction.
    f10 is used to weight the fractal in the z direction.
    i1 is the number of materials to use for the fractal distribution (defined according to the associated mixing model).
    This should be set to one if using a normal material instead of a mixing model.
    str1 is an identifier for the associated mixing model or material.
    str2 is an identifier for the fractal box itself.
    i2 is an optional parameter which controls the seeding of the random number generator used to create the fractals.
    By default (if you don’t specify this parameter) the random number generator will be seeded by trying to read data
    from /dev/urandom (or the Windows analogue) if available or from the clock otherwise.
    c1 is an optional parameter which can be y or n, used to switch on and off dielectric smoothing. If c1 is specified
    then a value for i2 must also be present.

    For example, to create an orthogonal parallelepiped with fractal distributed properties using a Peplinski mixing
    model for soil, with 50 different materials over a range of water volumetric fractions from 0.001 - 0.25, you
    should first define the mixing model using: #soil_peplinski: 0.5 0.5 2.0 2.66 0.001 0.25
    my_soil and then specify the fractal box using #fractal_box: 0 0 0 0.1 0.1 0.1 1.5 1 1 1 50 my_soil my_fractal_box.
    """
    txt = '#fractal_box: {} {} {} {} {} {} {} {} {} {} {} {} {}\n\n'.format(x0, y0, z0, x1, y1, z1,
                                                                           f7, f8, f9, f10, i1,
                                                                           str1, str2)
    return txt


def create_input_file(domain,
                      dx='0.002',
                      time_window='50e-9',
                      rx='0.002',
                      directoty=r"C:\Users\user\Documents\gprMax_project\dataset",
                      name_of_file="empty.in",
                      antenna_name='GSSI_400',
                      ):
    os.chdir(directoty)
    text = ''
    text += '#title: {}\n\n'.format(name_of_file[:-3], )
    text += '#domain: {} {} {}\n\n'.format(*domain)
    text += '#dx_dy_dz: {} {} {}\n\n'.format(dx, dx, dx)
    text += '#time_window: {}\n\n'.format(time_window)
    if antenna_name:
        if antenna_name == 'GSSI_400':
            text += '#python:\nfrom user_libs.antennas.GSSI import antenna_like_GSSI_400\n' \
                    'antenna_like_GSSI_400({} + current_model_run * {}, {}, {}, resolution={})\n' \
                    '#end_python:\n\n'.format('0.190',
                                              rx,
                                              round(domain[1] / 2, 3),
                                              round(domain[2] - 0.22, 3),
                                              dx)

    box_number = 0
    box_base = str(box_number) + '_soil'
    text += material_peplinski_box(material_name=box_base)
    text += geometry_peplinski_box(*domain, str1=box_base, str2=box_base+'_fractal')
    # print(text)
    with open(name_of_file, 'w') as f:
        f.write(text)

        # f.write('хуй')

import time
st=time.time()

# os.chdir(r"C:\Users\user\Documents\gprMax_project\dataset")
directoty = r"C:\Users\user\Documents\gprMax_project\dataset",
name_of_file="empty.in"
create_input_file(domain=[0.380, 0.380, 0.360])

api(r"C:\Users\user\Documents\gprMax_project\dataset\empty.in", geometry_only=True)

print("----%.2f----"%(time.time()-st))
