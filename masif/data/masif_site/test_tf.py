import tensorflow.compat.v1 as tf  
tf.compat.v1.disable_eager_execution()
tf.disable_v2_behavior()
if tf.test.gpu_device_name(): 

    print('Default GPU Device: {}'.format(tf.test.gpu_device_name()))

else:

   print("Please install GPU version of TF")
