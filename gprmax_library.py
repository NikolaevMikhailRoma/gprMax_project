import random
import os
from gprMax.gprMax import api
import h5py
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image
import scipy.ndimage as ndimage
from gprMax.receivers import Rx
from tools.plot_Ascan import mpl_plot as ascan_mpl_plot
from tools.outputfiles_merge import merge_files
import glob
from tools.plot_Bscan import get_output_data, mpl_plot
from scipy import interpolate


class GprMaxInputFile:
    def _init_dx(self, dx):
        if type(dx) != list:
            self.dx = [dx, dx, dx]
        else:
            self.dx = dx

    def _init_domain(self, domain):
        if len(domain) == 2:
            ### Здесь мы увеличиваем модель по горизонтали, чтобы вместить выносы антенны
            self.domain = [domain[0] + self.gpr_takeaway, domain[1] + self.gpr_takeaway, self.dx[2]]
        else:
            self.domain = domain

    def _init_random(self):
        rand = random.randint(0, 654123)
        return rand

    def __init__(self,
                 domain=[1.0, 1.0],
                 dx=0.002,
                 time_window='50e-9',
                 antenna_name='400',
                 wave_form='ricker',
                 rx_step=0.05,
                 directory=r"C:\Users\user\Documents\gprMax_project\dataset",
                 name_of_file="empty.in",
                 create_h5_geometry=False,
                 random_state=None,
                 random_pec_state=None,
                 availableoutputs=None

                 ):
        """
        :param domain[x, y, z]: (вверх, в глубину, в бок), выбираем площадь интересующей нас геометрии,
            дололнительно domain увеличивается в зависимости от антенны
        :param dx: разрешение модели по измерениям
        :param time_window: '50e-9' временная развертка, пока по умолчанию 50 ##
        :param rx: шаг измерения
        :param directory: с какой директорией мы работаем
        :param name_of_file: название файла импута
        :param antenna_name: название антенны или частота !!! в нынешнем варианте будем использовать 2D режим
        :param create_h5_geometry: создавать ли h5 геометрию
        :param random_state: random step
        """
        """GIPERPARAMS"""
        self.gpr_takeaway = 0.4

        self._init_dx(dx)
        self._init_domain(domain)

        self.time_window = time_window
        self.directory = directory
        self.name_of_file = name_of_file
        self.antenna_name = antenna_name
        self.wave_form = wave_form
        self.rx = rx_step
        self.num_of_steps = int(domain[0] / self.rx)
        self.create_h5_geometry = create_h5_geometry
        self.txt = None
        self.pic_size = 128
        self.availableoutputs = availableoutputs

        if random_state is None:
            self.random_state = self._init_random()
        if random_pec_state is None:
            random.seed(self._init_random())

    def material_peplinski_box(self, f1=0.5, f2=0.5, f3=2.0, f4=2.66,
                               f5=0.001, f6=0.25,
                               str1=None, ):
        """
        #soil_peplinski: f1 f2 f3 f4 f5 f6 str1

        :param f1: is the sand fraction of the soil
        :param f2: is the clay fraction of the soil.
        :param f3: is the bulk density of the soil in grams per centimetre cubed.
        :param f4: is the density of the sand particles in the soil in grams per centimetre cubed.
        :param f5: define a range for the volumetric water fraction of the soil (mix)
        :param f6: idem (max)
        :param str1: is an identifier for the soil
        :return: txt
        """

        if str1 is None:
            str1 = 'soil_peplinski_1'
        txt = f'#soil_peplinski: {f1} {f2} {f3} {f4} {f5} {f6} {str1}\n\n'

        if self.random_state:
            pass
        return txt

    def geometry_fractal_box(self, str1='soil_peplinski_1', str2=None, start_coordinate=None, end_coordinate=None,
                             f7=1.5, f8=3, f9=1, f10=1, i1=50, i2=False, c1=False):
        """
        Example
        #soil_peplinski: 0.5 0.5 2.0 2.66 0.001 0.25
        #fractal_box: 0 0 0 0.1 0.1 0.1 1.5 1 1 1 50 my_soil my_fractal_box.

        :param str1: name of peplinski material
        :param str2: name of box
        :param start_coordinate: first coordinate
        :param end_coordinate: end coordinate
        :param f7: is the fractal dimension which, for an orthogonal parallelepiped, should take values between zero and
                three.
        :param f8: is used to weight the fractal in the x direction.
        :param f9: is used to weight the fractal in the y direction.
        :param f10: is used to weight the fractal in the z direction.
        :param i1: is the number of materials to use for the fractal distribution (defined according to the associated
                mixing model)
        :param i2 is an optional parameter which controls the seeding of the random number generator used to create the
                fractals
        :param c1:is an optional parameter which can be y or n, used to switch on and off dielectric smoothing. If c1 is
                specified
        then a value for i2 must also be present.
        :return: text with box


        """
        if start_coordinate is None:
            x0 = y0 = z0 = 0
        else:
            x0, y0, z0 = start_coordinate[0], start_coordinate[1], start_coordinate[2]
        if end_coordinate is None:
            x1, y1, z1 = self.domain
        else:
            x1, y1, z1 = end_coordinate[0], end_coordinate[1], end_coordinate[2]

        if str2 is None:
            str2 = 'fractal_box_1'

        txt = f'#fractal_box: {x0} {y0} {z0} {x1} {round(y1, 2)} {z1} {f7} {f8} {f9} {f10} {i1} {str1} {str2}\n\n'
        if self.random_state:
            txt = f'#fractal_box: {x0} {y0} {z0} {x1} {y1} {z1} {f7} {f8} {f9} {f10} {i1} {str1} {str2} ' \
                  f'{self.random_state}\n\n'
        return txt

    def create_text(self):
        """
        сейчас сосредоточимся на геометрии
        :return:
        """

        self.txt = '#title: {}\n'.format(self.name_of_file[:-3])
        self.txt += '#domain: {} {} {}\n'.format(*self.domain)
        self.txt += '#dx_dy_dz: {} {} {}\n'.format(*self.dx)
        self.txt += f'#time_window: {self.time_window}\n'
        self.txt += f'#waveform: {self.wave_form} 1 {self.antenna_name}e6 my_ricker\n'
        self.txt += f'#hertzian_dipole: z {round(self.gpr_takeaway / 4, 2)} {round(self.domain[1] - self.gpr_takeaway, 2)} 0 my_ricker\n'
        self.txt += f'#rx: {round(self.gpr_takeaway * 3 / 4, 2)} {round(self.domain[1] - self.gpr_takeaway, 2)} 0\n'
        self.txt += f'#src_steps: {self.rx} 0 0\n#rx_steps: {self.rx} 0 0\n\n'
        # selt.txt += '#num_threads: 8'
        # todo убрал блок пеплинского
        self.txt += self.material_peplinski_box()
        self.txt += self.geometry_fractal_box(end_coordinate=[self.domain[0], self.domain[1] - self.gpr_takeaway,
                                                              self.domain[2]],
                                              )

        # # self.txt += '#material: 6 0 1 0 half_space\n'
        # # self.txt += '#box: 0 0 0 {} {} {} half_space\n\n'.format(self.domain[0],
        # #                                                          round_custom(self.domain[1] - self.gpr_takeaway, self.dx[0]),
        #                                                          self.domain[2])


        # # self.create_h5_geometry = True
        # if self.create_h5_geometry:
        #     self.txt += f'#geometry_objects_write: {round(self.gpr_takeaway / 2, 2)} 0 0 ' \
        #                 f'{round(self.domain[0] - self.gpr_takeaway / 2, 2)} {round(self.domain[1] - self.gpr_takeaway, 2)} ' \
        #                 f'{self.domain[2]} geom'

        for i in range(random.randint(1, 4)):
            self.txt += self.add_pec()

        ### Check orientation
        # self.txt += f'#cylinder: {self.gpr_takeaway/2} {0.1} {0} ' \
        #             f'{self.gpr_takeaway/2} {0.1} {self.domain[2]} {0.05} pec\n'
        # self.txt += f'#cylinder: {self.domain[0]/2} {0.9} {0} {self.domain[0]/2} {0.9} {self.domain[2]} {0.05} pec\n'

    def add_pec(self):
        x = round_custom(random.uniform(self.gpr_takeaway*3/4, self.domain[0]-self.gpr_takeaway*3/4), self.dx[0])
        y = round_custom(random.uniform(0, self.domain[1] - self.gpr_takeaway - 0.2), self.dx[0])
        radius = round_custom(random.uniform(self.dx[0]*5, self.dx[0]*10), self.dx[0])

        text_pec = f'#cylinder: {x} {y} {0} {x} {y} {self.domain[2]} {radius} pec\n'
        return text_pec

    def save_input_file(self):
        os.chdir(self.directory)
        with open(self.name_of_file, 'w') as f:
            f.write(self.txt)

    def start_modelling(self, n=None):
        if n is None:
            api(self.directory + '\\' + self.name_of_file,
                n=int(round((asdf.domain[0] - asdf.gpr_takeaway) / asdf.rx, 0)), gpu={0}, geometry_only=False)
        else:
            api(self.directory + '\\' + self.name_of_file, n=n, gpu={0}, geometry_only=False)

    def create_geometry_file(self):
        os.chdir(self.directory)
        name_of_geometry_file = self.name_of_file[:-3] + '_geometry' + self.name_of_file[-3:]
        # add mini space for adecvation geometry materials file
        text_for_geometry = self.txt + f'#geometry_objects_write: {round(self.gpr_takeaway / 2, 2)} 0 0 ' \
                                       f'{round(self.domain[0] - self.gpr_takeaway / 2, 2)} ' \
                                       f'{round(self.domain[1] - self.gpr_takeaway + self.dx[2], 3)} ' \
                                       f'{self.domain[2]} geom'

        with open(name_of_geometry_file, 'w') as f:
            f.write(text_for_geometry)

        api(name_of_geometry_file, geometry_only=True)

    def geometry_file_manipulation(self, save_to=None, show=False):
        with open(self.directory + '\geom_materials.txt', 'r') as f:
            material_file = f.readlines()
        material_file = ''.join([str(elem) for elem in material_file]).split('#material: ')[1:]

        relative_permittivity_materials_file = np.array([float(x.split(' ')[0]) for x in material_file])
        replacement_perm = lambda t: relative_permittivity_materials_file[t]
        vfunc_perm = np.vectorize(replacement_perm)

        conductivity_materia_file = np.array([float(x.split(' ')[1]) for x in material_file])
        conductivity_materia_file[0] = 1024
        replacement_cond = lambda t: conductivity_materia_file[t]
        vfunc_cond = np.vectorize(replacement_cond)

        f1 = h5py.File(self.directory + "\geom.h5", 'r')
        f1 = np.array(f1['data'][:, :, :])

        f1 = ndimage.rotate(f1, 90, reshape=True)

        geom_perm = vfunc_perm(f1)
        geom_cond = vfunc_cond(f1)
        # print(np.unique(geom_cond))

        geom_file = np.stack((geom_perm, geom_cond), axis=-1)
        geom_file = geom_file.reshape([geom_file.shape[0], geom_file.shape[1], 2])

        res = geom_file[:, :, :]
        plt.imshow(res[:, :, 0], interpolation='nearest')

        print(res.shape)
        if save_to == 'png':
            plt.savefig(self.directory + '\\' + self.name_of_file[:-3] + '.png')

        elif save_to == 'np':
            np.save(self.directory + '\\' + 'geom_' + self.name_of_file[:-3], res)
        if show:
            plt.show()
        # plt.colorbar()
        plt.close()
        return res

    def snow_np_pic(self):
        npy = np.load(self.directory + '\\'+'geom_' + self.name_of_file[:-3] + '.npy')
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(npy[:, :, 0])
        ax[1].imshow(npy[:, :, 1])
        plt.show()
        plt.close()

    def ascan_plot(self, filename):
        filename = self.directory + '\\' + filename
        # api(filename, n=1, geometry_only=False)
        outputs = Rx.defaultoutputs
        plt = ascan_mpl_plot(filename, outputs, fft=False)
        plt.show()
        plt.close()

    def merge_data(self, removefiles=False, version=2):
        basefilename = self.directory + '\\' + self.name_of_file[:-3]
        if version == 1:
            merge_files(basefilename=basefilename, removefiles=removefiles)
        elif version == 2:
            files = glob.glob(basefilename + '[0-9]*.out')
            outputfiles = [filename for filename in files if '_merged' not in filename]
            modelruns = len(outputfiles)
            # print(outputfiles)
            outputfile = basefilename + '_merged.out'
            # Combined output file
            fout = h5py.File(outputfile, 'w')
            for model in range(modelruns):
                fin = h5py.File(basefilename + str(model + 1) + '.out', 'r')
                nrx = fin.attrs['nrx']

                # Write properties for merged file on first iteration
                if model == 0:
                    fout.attrs['Title'] = fin.attrs['Title']
                    fout.attrs['gprMax'] = '32123'
                    fout.attrs['Iterations'] = fin.attrs['Iterations']
                    fout.attrs['dt'] = fin.attrs['dt']
                    fout.attrs['nrx'] = fin.attrs['nrx']
                    for rx in range(1, nrx + 1):
                        path = '/rxs/rx' + str(rx)
                        grp = fout.create_group(path)
                        if self.availableoutputs:
                            availableoutputs = self.availableoutputs
                        else:
                            availableoutputs = list(fin[path].keys())
                        for output in availableoutputs:
                            grp.create_dataset(output, (fout.attrs['Iterations'], modelruns),
                                               dtype=fin[path + '/' + output].dtype)

                # For all receivers
                for rx in range(1, nrx + 1):
                    path = '/rxs/rx' + str(rx) + '/'
                    if self.availableoutputs:
                        availableoutputs = self.availableoutputs
                    else:
                        availableoutputs = list(fin[path].keys())
                    # For all receiver outputs
                    for output in availableoutputs:
                        fout[path + '/' + output][:, model] = fin[path + '/' + output][:]

                fin.close()

            fout.close()
            if removefiles:
                for model in range(modelruns):
                    file = basefilename + str(model + 1) + '.out'
                    os.remove(file)

    def default_bscan_plot(self, output=None):
        basefilename = self.directory + '\\' + self.name_of_file[:-3]
        filename = basefilename + '_merged.out'
        rxnumber = 1
        if output is None:
            rxcomponent = 'Ez'
        else:
            rxcomponent = self.availableoutputs

        outputdata, dt = get_output_data(filename, rxnumber, rxcomponent)
        plt = mpl_plot(outputdata=outputdata, dt=dt, rxnumber=rxnumber, rxcomponent=rxcomponent, filename=filename)
        plt.show()
        plt.close()

    def create_bsan_npy(self, removefiles=False, len_of_trace=False):
        basefilename = self.directory + '\\' + self.name_of_file[:-3]

        files = glob.glob(basefilename + '[0-9]*.out')
        outputfiles = [filename for filename in files if '_merged' not in filename]
        modelruns = len(outputfiles)

        # Combined output file
        for model in range(modelruns):
            fin = h5py.File(basefilename + str(model + 1) + '.out', 'r')
            nrx = fin.attrs['nrx']

            # Write properties for merged file on first iteration
            if model == 0:
                path = '/rxs/rx1/' + self.availableoutputs[0]

                npy_out = np.array(fin[path][:])
                npy_out = npy_out.reshape([1, len(npy_out)])
                # print(npy_out.shape)
            else:
                npy_out = np.append(npy_out, [fin[path][:]], axis=0)
            fin.close()

        if len_of_trace:
            shape_of_array = npy_out.shape
            step = shape_of_array[1] / len_of_trace
            for i in range(shape_of_array[0]):

                y = npy_out[i]
                x = range(0, shape_of_array[1])
                f = interpolate.interp1d(x, y)
                xnew = np.arange(0, shape_of_array[1], step)
                ynew = f(xnew)
                if i == 0:
                    npy_out_new = np.array(ynew).reshape([1, len_of_trace])
                else:
                    npy_out_new = np.append(npy_out_new, [ynew], axis=0)
            npy_out = npy_out_new

        print(npy_out.shape)
        np.save(self.directory + '\\' + 'bscan_' + self.name_of_file[:-3], npy_out)
        # print(npy_out.shape)
        if removefiles:
            for model in range(modelruns):
                file = basefilename + str(model + 1) + '.out'
                os.remove(file)


def round_custom(num, step):
    return round(round(num/step) * step, 2)


if __name__ == '__main__':
    # random.seed() # параметр отвечает за рандом пеков
    asdf = GprMaxInputFile(rx_step=0.008, antenna_name='800', time_window='30e-9', availableoutputs=['Ez'], dx=0.008)
    asdf.create_text()
    asdf.save_input_file()
    asdf.create_geometry_file()
    asdf.geometry_file_manipulation(save_to='np', show=True)
    asdf.geometry_file_manipulation(save_to='png', show=False)

    asdf.snow_np_pic()
    asdf.start_modelling(n=128)
    asdf.merge_data(removefiles=False, )
    asdf.default_bscan_plot()
    asdf.create_bsan_npy(len_of_trace=512, removefiles=True)

    # asdf.save_input_file()
    # asdf.merge_data(removefiles=False)
    # asdf.start_modelling()
    # asdf.create_geometry_file()
    # asdf.show_geometry_file(save_to=None)
    # asdf.ascan_plot('empty14.out')
