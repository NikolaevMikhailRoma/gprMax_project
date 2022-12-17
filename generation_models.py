import gprmax_library

if __name__ == '__main__':
    for i in range(1000, 2000):
        empty_name = "%06d" % (i)
        asdf = gprmax_library.GprMaxInputFile(directory=r"C:\Users\user\Documents\gprMax_project\dataset_2",
                                              name_of_file=empty_name+'.in',
                                              rx_step=0.008, antenna_name='800', time_window='30e-9',
                                              availableoutputs=['Ez'], dx=0.008)
        asdf.create_text()
        asdf.save_input_file()
        asdf.create_geometry_file()
        asdf.geometry_file_manipulation(save_to='np', show=False)
        # asdf.snow_np_pic()
        try:
            asdf.start_modelling(n=128)
        except:
            asdf.start_modelling(n=128)

        # try:
        #
        # except:
        #     asdf.start_modelling(n=128)

        # asdf.merge_data(removefiles=False,)
        asdf.create_bsan_npy(len_of_trace=512, removefiles=True)

