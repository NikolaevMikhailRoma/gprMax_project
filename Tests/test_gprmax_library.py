from gprmax_library import GprMaxInputFile
import pytest
import unittest


class TestInputFile2D(unittest.TestCase):
    def setUp(self, domain=[1.00, 1.00, 0.002], dx=[0.002, 0.002, 0.002], time_window='50e-9',
              antenna_name='400', wave_form='ricker', rx_step='0.002',
              directory=r"C:\Users\user\Documents\gprMax_project\dataset", name_of_file="empty.in",
              create_h5_geometry=False, random_state=False):
        self.dx = dx
        self.domain = domain
        self.time_window = time_window
        self.directory = directory
        self.name_of_file = name_of_file
        self.antenna_name = antenna_name
        self.wave_form = wave_form
        self.rx_step = rx_step
        self.create_h5_geometry = create_h5_geometry
        self.random_state = random_state
        # self.input_file = GprMaxInputFile(domain=domain, dx=dx, time_window=time_window, antenna_name=antenna_name,
        #                                   wave_form=wave_form, rx=rx, directory=directory,
        #                                   name_of_file=name_of_file,
        #                                   create_h5_geometry=create_h5_geometry)

    def test_gpr_max_input_file(self):
        print("\n@@@@@@@@Test gpr_max_input_file@@@@@@@@@")
        input_file = GprMaxInputFile()
        self.assertEqual(input_file.domain, self.domain)
        assert input_file.domain == self.domain
        assert input_file.dx == self.dx
        assert input_file.time_window == self.time_window
        assert input_file.antenna_name == self.antenna_name
        assert input_file.wave_form == self.wave_form
        assert input_file.rx_step == self.rx_step
        assert input_file.directory == self.directory
        assert input_file.name_of_file == self.name_of_file
        assert input_file.create_h5_geometry == self.random_state

    def test_material_peplinski_box(self):
        print("\n@@@@@@@@test_material_peplinski_box@@@@@@@@@")
        input_file = GprMaxInputFile()
        txt = input_file.material_peplinski_box()
        self.assertEqual(txt, f'#soil_peplinski: 0.5 0.5 2.0 2.66 0.001 0.25 1\n\n')

    def test_fractal_box(self):
        print("\n@@@@@@@@ test_fractal_box @@@@@@@@@")

        input_file = GprMaxInputFile()
        txt = input_file.geometry_fractal_box()
        self.assertEqual(txt, '#fractal_box: 0 0 0 1.0 1.0 0.002 1.5 1 1 1 50 1 fractal_box_1\n\n')



#
#
# def test_geometry_peplinski_box():
#     assert False
