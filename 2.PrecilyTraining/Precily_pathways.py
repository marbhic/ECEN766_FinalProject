##Importing libraries1923
import numpy as np
import pandas as pd
import keras_tuner as kt
import tensorflow as tf
from tensorflow import keras

tf.config.threading.set_inter_op_parallelism_threads(1)
tf.config.threading.set_intra_op_parallelism_threads(1)
np.random.seed(2)

##loading 90% Training data
train_set = pd.read_csv('C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\training_data_attempt2.csv')
directory = 'C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\NewDataTest'
project_name = 'Pathways'


##Definining hyper parameters 
layers_range = (2, 6)
units_range = (128, 256, 4)
lr_values = [1e-3,1e-4,1e-5]

##Define model
def model_builder(hp):
  model = keras.Sequential()
  model.add(keras.layers.Dense(1923, input_dim = 1923, activation = 'relu'))
  model.add(keras.layers.Dense(units = 512, activation = 'relu'))
  model.add(
            keras.layers.Dropout(
                hp.Float(
                    'dropout',
                    min_value=0.1,
                    max_value=0.5,
                    default=0.1,
                    step=0.1)
            )
        )

  for i in range(hp.Int('layers', layers_range[0], layers_range[1])):
    model.add(keras.layers.Dense(units=hp.Int('units_' + str(i),  
                                min_value=units_range[0], max_value=units_range[1], 
                                step=units_range[2]), activation='relu'))
    model.add(
            keras.layers.Dropout(
                hp.Float(
                    'dropout',
                    min_value=0.1,
                    max_value=0.5,
                    default=0.1,
                    step=0.1)
            )
        )

    
  model.add(keras.layers.Dense(1))
  hp_learning_rate = hp.Choice('learning_rate', values=lr_values)
  model.compile(optimizer=keras.optimizers.Adam(learning_rate=hp_learning_rate),
                loss='mean_squared_error') 
  return model

##Perfrom cell line wise split
#def cell_line_split(data, k):
  #CL = data['CELL_LINE_NAME'].unique()
  #np.random.shuffle(CL)
  #A = set(CL[:int(k*len(CL))])
  #B = set(CL[int(k*len(CL)):])
  #train, test = [], []
  
  #for i in df.to_numpy():
      #if i[0] in A:
          #train.append(i)
      #else:
          #test.append(i)
  #train = pd.DataFrame(train)
  test = pd.DataFrame(test)
  #return train, test

##If complete data is being used use perform cell line wise split using cell_line_split function
#train_set, test_set = cell_line_split(df, 0.9)
#df = None
#train_set.to_csv("Train_data.csv",index=False)
#test_set.to_csv("Test_data.csv",index=False)
#test_set = None

CL_x = train_set[train_set.columns[0]].unique()
CL_x = list(CL_x)

##Perform hyperparameter tuning
for i in range(5):
    A = set( CL_x[:i*len(CL_x)//5] + CL_x[(i+1)*len(CL_x)//5:] )
    B = set( CL_x[i*len(CL_x)//5:(i+1)*len(CL_x)//5] )
    
    train, val = [], []
    for j in train_set.to_numpy():
        if j[0] in A: 
            train.append(j)
        else:
            val.append(j)

    train = pd.DataFrame(train)
    val = pd.DataFrame(val)
    # train.to_csv(path+"Train_Set_"+str(i+1)+".csv",index=False)
    # val.to_csv(path+"Validation_Set_"+str(i+1)+".csv",index=False)
    X_train = train.iloc[: , 2:-1]
    Y_train = train.iloc[: , -1:]
    X_val = val.iloc[: , 2:-1]
    Y_val = val.iloc[: , -1:]
    train, test = None, None
    X_train.to_csv("Train_Set_"+str(i+1)+".csv",index=False)
    Y_train.to_csv("Train_Set_truth_"+str(i+1)+".csv",index=False)
    X_val.to_csv("Val_Set_"+str(i+1)+".csv",index=False)
    Y_val.to_csv("Val_Set_truth_"+str(i+1)+".csv",index=False)
    
    tuner = kt.Hyperband(model_builder, # the hypermodel
                    objective='val_loss', # objective to optimize
                    max_epochs=30,
                    factor=3, 
                    directory=directory, # directory to save logs 
                    project_name=project_name+str(i+1))
    
    stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)
    tuner.search(X_train, Y_train, epochs=30, validation_data = (X_val, Y_val), callbacks=[stop_early])
    best_hp=tuner.get_best_hyperparameters()[0]
    best_model = tuner.get_best_models()[0]
    
    # Build the model with the optimal hyperparameters
    h_model = tuner.hypermodel.build(best_hp)
    h_model.fit(X_train, Y_train, epochs=50, verbose = 1, batch_size = 128, validation_data = (X_val, Y_val))
    h_model.save('C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\precily_cv_'+str(i+1)+'.hdf5')
    h_model.save('C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\precily_cv_'+str(i+1)+'.keras')
    #h_model = None

    # Train using whole model
    train2 = pd.read_csv("C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\training_data_attempt2.csv")
    X_train1 = train2.iloc[: , 2:-1]
    Y_train1 = train2.iloc[: , -1:]
    h_model.fit(X_train1, Y_train1, verbose = 1, epochs=50, batch_size = 128)
    h_model.save('C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\Model'+str(i+1)+'.hdf5')
    h_model.save('C:\\Users\\danbeltran\\Documents\\ECEN766\\20240419PrecilyTraining\\Model'+str(i+1)+'.keras')
    
    

##############################Training of complete dataset using hyper-tuned models###################################
######################################################################################################################

# ###Load training data
# train = pd.read_csv("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\Fig1\\Fig1c\\Fig1c_Precily_pathways\\training_data_attempt2.csv")
# X_train = train.iloc[: , 2:-1]
# Y_train = train.iloc[: , -1:]

# ###Load hyper-tuned models obtained from the previous step and fit training data
# for i in range(5):
#     h_model = tf.keras.models.load_model("C:\\Users\\marbhic\\Downloads\\Precily-main\\Precily-main\\"+'precily_cv_'+str(i+1)+'.hdf5')
#     print(h_model)
#     print("hi")
#     h_model.fit(X_train, Y_train, verbose = 1, epochs=50, batch_size = 128)
#     h_model.save('Model'+str(i+1)+'.hdf5')
#     h_model = None