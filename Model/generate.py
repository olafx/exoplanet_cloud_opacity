'''
Generate, train, and store models along with their performance criteria.
'''

import pytz

#   SETTINGS

#   File for data.
filename_data          = '../Data/data.h5'
data_main_group_name   = 'main'
#   Which datasets to use in training and testing.
data_names             = [1, 2]
#   Ratio of training data to testing data.
train_to_test_ratio    = 0.95
#   Number of epochs.
n_epochs               = 10
#   File for model performance criteria.
filename_models        = '../Data/models.h5'
models_main_group_name = 'models'
#   Folder where model will be stored.
foldername_models      = '../Data/Models'
#   Computer timezone. (Needed to keep consistent timestamps.)
source_zone            = pytz.timezone('Europe/Brussels')
#   Store model or not?
store                  = False

import tensorflow as tf
import numpy as np
import h5py
from datetime import datetime
import math

target_zone = pytz.timezone('UTC')



#   Abstractify obtaining and preprocessing testing and training data.
class Data:

    def __init__(self):
        #   Opening filename_data and verifying that data_main_group_name exists.
        if data_main_group_name not in (fp := h5py.File(filename_data, 'r')).keys():
            raise Exception(f"group '{data_main_group_name}' does not exist in file, or file does not exist")
        #   Extracting data requested in data_names.
        group, data_in, data_out = fp[data_main_group_name], [], []
        for name in data_names:
            data_in.append(group[str(name)]['in']), data_out.append(group[str(name)]['out'])
        data_in, data_out = np.concatenate(data_in), np.concatenate(data_out)
        fp.close()
        #   Mixing data. This is necessary since various data_names may have differing settings.
        np.random.seed(None) # Runtime seeded by entropy source.
        state_initial = np.random.get_state()
        np.random.shuffle(data_in)
        np.random.set_state(state_initial)
        np.random.shuffle(data_out)
        self.input, self.output = data_in, data_out
        print(f'initialized {len(self)} samples')

    def __len__(self):
        #   Number of data points.
        return self.input.shape[0]

    def __preprocess_input(self):
        #   Reducing dimensionality of mixture ratios to achieve linear independence.
        self.input = self.input[:, : -1]
        #   Wavelengths and dust sizes are scaled logarithmatically, followed by theoretical minmax
        #   normalization. (Based on theoretical min and max from RNG settings as opposed to actual
        #   min and max in data.)
        #   Dust size sd is also minmax scaled, but not treated logarithmically. (And has implicit 0
        #   min.)
        #   Different data_names may have different settings, so finding correct min and max is
        #   not trivial.
        fp = h5py.File(filename_data, 'r')
        #   Finding min and max for all properties that need scaling.
        range_wavelength = [math.inf, -math.inf]
        range_size       = [math.inf, -math.inf]
        max_size_sd      = -math.inf
        for name in data_names:
            group = fp[data_main_group_name + '/' + str(name)]
            if ((min_wavelength := group.attrs['wavelength minimum']) < range_wavelength[0]):
                range_wavelength[0] = min_wavelength
            if ((max_wavelength := group.attrs['wavelength maximum']) > range_wavelength[1]):
                range_wavelength[1] = max_wavelength
            if ((min_size := group.attrs['dust size minimum']) < range_size[0]):
                range_size[0] = min_size
            if ((max_size := group.attrs['dust size maximum']) > range_size[1]):
                range_size[1] = max_size
            if ((max_size_sd := group.attrs['dust size standard deviation maximum']) > max_size_sd):
                max_size_sd = max_size_sd
        fp.close()
        #   Min max and logarithmic scaling.
        self.input[:, 0] = (np.log(self.input[:, 0])    - np.log(range_wavelength[0])) /\
                           (np.log(range_wavelength[1]) - np.log(range_wavelength[0]))
        self.input[:, 1] = (np.log(self.input[:, 1])    - np.log(range_size[0])) /\
                           (np.log(range_size[1])       - np.log(range_size[0]))
        self.input[:, 2] = self.input[:, 2] / max_size_sd
        print('minmax scaler settings')
        print('\trange wavelength', range_wavelength)
        print('\trange size', range_size)
        print('\tmaximum size sd', max_size_sd)
        print('preprocessed input')
        return np.array(range_wavelength), np.array(range_size), max_size_sd

    def __preprocess_output(self):
        #   No out preprocessing is done.
        print('preprocessed output')

    def preprocess(self):
        self.__preprocess_output()
        return self.__preprocess_input()

    def train(self):
        #   Get training data, split in physical input properties, chemical input properties, and
        #   output opacity properties.
        self.__train_to_test_ratio_validity()
        end = int(len(self) * train_to_test_ratio)
        return self.input[: end, : 3], self.input[: end, 3 :], self.output[: end]

    def test(self):
        #   Get testingdata, split in physical input properties, chemical input properties, and
        #   output opacity properties.
        self.__train_to_test_ratio_validity()
        start = int(len(self) * train_to_test_ratio)
        return self.input[start :, : 3], self.input[start :, 3 :], self.output[start :]

    def __train_to_test_ratio_validity(self):
        if train_to_test_ratio <= 0 or train_to_test_ratio >= 1:
            raise Exception('train_to_test_data not valid')



#   Abstractify storage for model and model properties, including performance criteria, training
#   datasets, and information required to reverse the preprocessing.
class Models:

    def __init__(self):
        #   Opening filename_models and verifying that models_main_group_name exists.
        self.fp = h5py.File(filename_models, 'a')
        if models_main_group_name not in self.fp.keys():
            raise Exception(f"group '{models_main_group_name}' does not exist in file, or file does not exist")
        self.main_group = self.fp[models_main_group_name]
        #   Verifying consistency; there should exist a previous model, and the new model shouldn't
        #   override another.
        if 0 < (n_models := len(self)):
            if str(n_models) not in self.main_group:
                raise Exception(f"'{n_models}' not in '{models_main_group_name}'")
            if str(n_models + 1) in self.main_group:
                raise Exception(f"'{n_models + 1}' already in '{models_main_group_name}'")

    def __del__(self):
        self.fp.close()

    def __len__(self):
        #   Number of models.
        return len(self.main_group)

    #   Adding a model along with its performance during training and testing.
    #   Also need to add minmax properties to know how to undo minmax scaling.
    def add(self, model, metrics, metrics_train, metrics_test, range_wavelength, range_size, max_size_sd):
        #   Add model to models file and set timestamp.
        (group := self.main_group.create_group(name := str(len(self) + 1))).attrs['creation date'] =\
            str(source_zone.localize(datetime.now()).astimezone(target_zone))
        #   Add performance metrics.
        group['metric names']           = metrics
        group['metrics train']          = metrics_train
        group['metrics test']           = metrics_test
        group['training file']          = filename_data
        group['training sets']          = data_names
        group.attrs['range_wavelength'] = range_wavelength
        group.attrs['range_size']       = range_size
        group.attrs['max_size_sd']      = max_size_sd
        #   Storing the actual model object in foldername_models.
        model.save(foldername_models + '/' + name + '.h5')



#   Abstractify model (network).
class Model:

    def __init__(self):
        #   Split of physical and chemical data. (Allows for more flexibility.)
        input0 = tf.keras.Input(shape=(3,))
        input1 = tf.keras.Input(shape=(15,))
        x = tf.keras.Model(input0, input0)
        y = tf.keras.Model(input1, input1)
        z = tf.keras.layers.concatenate([x.output, y.output])
        z = tf.keras.layers.Dense(12, activation='gelu')(z)
        z = tf.keras.layers.LayerNormalization()(z)
        z = tf.keras.layers.Dense(9, activation='gelu')(z)
        z = tf.keras.layers.Dense(8, activation='gelu')(z)
        z = tf.keras.layers.Dense(7, activation='gelu')(z)
        z = tf.keras.layers.Dense(6, activation='gelu')(z)
        z = tf.keras.layers.Dense(5, activation='linear')(z)
        z = tf.keras.layers.Dense(2, activation='linear')(z)
        self.model = tf.keras.Model([x.input, y.input], z)

    def compile(self):
        self.model.compile(
            optimizer='nadam',
            loss=tf.keras.losses.MeanAbsoluteError(),
            metrics=[tf.metrics.MeanAbsoluteError(), tf.metrics.LogCoshError(), tf.metrics.MeanSquaredError()]
        )



def main():

    #   Set up data.
    data = Data()
    data.train_to_test_ratio = train_to_test_ratio
    range_wavelength, range_size, max_size_sd = data.preprocess()

    #   Set up models.
    models = Models()

    #   Set up model.
    model = Model()
    model.compile()
    #   Only need the actual TF model from here on.
    model = model.model

    #   Acquiring training and testing data.
    (train_in_physical, train_in_chemical, train_out), (test_in_physical, test_in_chemical, test_out) = data.train(), data.test()

    #   Training and testing, and saving results.
    results_train, results_test = [], []
    for i in range(n_epochs):
        results_train.append(model.fit([train_in_physical, train_in_chemical], train_out))
        results_test.append(model.evaluate([test_in_physical, test_in_chemical], test_out))

    #   Getting metric names. Fixing wrong string type.
    metrics = np.array(list(results_train[0].history), dtype=np.str).astype('S')
    #   Getting values of training metrics. Original data structure is strange, requires fixing.
    
    metrics_train = [[result[0] for result in list(result.history.values())] for result in results_train]
    
    #   Getting values of testing metrics.
    metrics_test = list(results_test)

    #   Store this model and performance history.
    if store:
        models.add(model, metrics, metrics_train, metrics_test, range_wavelength, range_size, max_size_sd)

if __name__ == '__main__':
    main()
