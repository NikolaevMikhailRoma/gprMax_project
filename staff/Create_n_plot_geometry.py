from gprMax.gprMax import api
import h5py
import numpy as np

# filename = 'user_models\cylinder_Bscan_GSSI_1500.in'
# filename_out = filename[:-3]+'.out'
# if not os.path.exists(filename_out):
#     api(filename, n=1, geometry_only=False, gpu={0})
#
# outputs = Rx.defaultoutputs
# outputs = ['Ez']
# plt = mpl_plot(filename_out, outputs, fft=False)
# plt.show()
from staff.Create_input_file import create_input_file


if __name__ == '__main__':

    #### создаем инпут файлик и геометрию
    directoty=r"C:\Users\user\Documents\gprMax_project\dataset"
    name_of_file="empty.in"

    input_text = create_input_file(domain=[1, 0.380, 1],
                      rx='0.002',
                      antenna_name='GSSI_400',
                      directoty=directoty,
                      name_of_file=name_of_file,
                      create_h5_geometry=True,
                      write_file=False)
    with open(directoty+'\\'+name_of_file, 'r') as f:
        lines = f.readlines()

    if ''.join(lines) == input_text:
            print('Акуальный файл')
    else:
        print('Create_geometry_file')
    api(r"C:\Users\user\Documents\gprMax_project\dataset\empty.in", geometry_only=True)
    #### первая часть сделанна

    #### разбираемся с h5 файлом

    f1 = h5py.File(directoty+'\\'+'geom.h5', 'r')
    print(list(f1.keys()))

    # считываем материалы
    with open(r"/dataset/geom_materials.txt", 'r') as f:
        material_file = f.readlines()
    material_file = material_file[0:-1:2]
    material_file = np.array([float(x.split(' ')[1]) for x in material_file])
    # print(material_file)


    # пробуем считать геометрию
    f1 = f1['data'][:,95,:]
    # print(f1)
    # print(np.unique(f1[:,:,:]))



    replacement = lambda t: material_file[t]
    vfunc = np.vectorize(replacement)
    f1 = vfunc(f1)

    print(f1)
    print('конец')

    from matplotlib import pyplot as plt

    plt.imshow(f1, interpolation='nearest')
    plt.colorbar()
    plt.show()