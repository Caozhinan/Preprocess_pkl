#!/usr/bin/env python3  
# -*- coding: utf-8 -*-  
  
import sys  
import os  
import time  
import numpy as np  
import argparse  
import importlib  
from default_config.masif_opts import masif_opts  
  
# Apply mask to input_feat  
def mask_input_feat(input_feat, mask):  
    mymask = np.where(np.array(mask) == 0.0)[0]  
    return np.delete(input_feat, mymask, axis=2)  
  
def load_hdf5_weights(learning_obj, hdf5_weight_path):  
    """加载HDF5权重文件到TensorFlow模型"""  
    import h5py  
      
    print("使用HDF5权重文件: {}".format(hdf5_weight_path))  
      
    with h5py.File(hdf5_weight_path, 'r') as f:  
        # 获取所有可训练变量  
        trainable_vars = learning_obj.session.graph.get_collection('trainable_variables')  
          
        # 打印模型参数信息（类似你日志中的输出）  
        # total_params = 0  
        # for var in trainable_vars:  
        #     print(var)  
        #     param_count = np.prod(var.get_shape().as_list())  
        #     print(param_count)  
        #     total_params += param_count  
          
        # print("Total number parameters: {}".format(total_params))  
          
        # 从HDF5文件加载权重  
        for var in trainable_vars:  
            var_name = var.name.replace(':0', '')  
            if var_name in f:  
                weight_data = f[var_name][:]  
                learning_obj.session.run(var.assign(weight_data))  
                # print("已加载权重: {}".format(var_name))  
            else:  
                print("警告: 在HDF5文件中未找到权重: {}".format(var_name))  
  
def main():  
    parser = argparse.ArgumentParser(description='Batch fingerprint generation for multiple complexes with HDF5 weights')  
    parser.add_argument('--precomputed_dirs', nargs='+', required=True,   
                       help='List of precomputed directories')  
    parser.add_argument('--output_dirs', nargs='+', required=True,  
                       help='List of output directories')  
    parser.add_argument('--ppi_pair_ids', nargs='+', required=True,  
                       help='List of PPI pair IDs')  
    parser.add_argument('--custom_params_file', required=False,  
                       help='Custom parameters file path')  
      
    args = parser.parse_args()  
      
    # 验证输入参数长度一致  
    if not (len(args.precomputed_dirs) == len(args.output_dirs) == len(args.ppi_pair_ids)):  
        print("错误: precomputed_dirs, output_dirs, 和 ppi_pair_ids 的长度必须相同")  
        sys.exit(1)  
      
    # 加载参数配置  
    params = masif_opts["ppi_search"]  
      
    # 尝试加载自定义参数  
    if args.custom_params_file and os.path.exists(args.custom_params_file):  
        try:  
            custom_params = importlib.import_module(args.custom_params_file.replace('.py', '').replace('/', '.'), package=None)  
            custom_params = custom_params.custom_params  
            for key in custom_params:  
                print("设置 {} 为 {}".format(key, custom_params[key]))  
                params[key] = custom_params[key]  
        except Exception as e:  
            # print("警告: 无法加载自定义参数文件 {}: {}".format(args.custom_params_file, str(e)))  
            print("将使用默认参数。")  
      
    # 一次性加载预训练的神经网络模型  
    print("加载预训练的神经网络模型...")  
    from masif_modules.MaSIF_ppi_search import MaSIF_ppi_search  
      
    learning_obj = MaSIF_ppi_search(  
        params["max_distance"],  
        n_thetas=16,  
        n_rhos=5,  
        n_rotations=16,  
        idx_gpu="/gpu:0",  
        feat_mask=params["feat_mask"],  
    )  
      
    # 使用HDF5权重文件  
    hdf5_weight_dir = "/es01/paratera/sce0413/czn/preprocess_pkl/masif/data/masif_ppi_search/nn_models/sc05/all_feat/weight"    
    hdf5_weight_path = os.path.join(hdf5_weight_dir, "weights_12A_0129.hdf5")  
      
    if os.path.exists(hdf5_weight_path):  
        load_hdf5_weights(learning_obj, hdf5_weight_path)  
        print("神经网络模型加载成功。")  
    else:  
        print("错误: HDF5权重文件不存在: {}".format(hdf5_weight_path))  
        sys.exit(1)  
      
    # 导入计算描述符的函数  
    from masif_modules.train_ppi_search import compute_val_test_desc  
      
    print("开始批量处理 {} 个复合物...".format(len(args.precomputed_dirs)))  
      
    # 循环处理所有复合物  
    for count, (precomputed_dir, output_dir, ppi_pair_id) in enumerate(zip(args.precomputed_dirs, args.output_dirs, args.ppi_pair_ids)):  
        print("处理复合物 {}/{}: {}".format(count + 1, len(args.precomputed_dirs), ppi_pair_id))  
          
        # 检查预计算目录是否存在  
        if not os.path.exists(precomputed_dir):  
            print("警告: 预计算目录不存在: {}".format(precomputed_dir))  
            continue  
          
        # print("从目录加载预计算特征: {}".format(precomputed_dir))  
          
        # 创建输出目录  
        out_desc_dir = os.path.join(output_dir, "descriptors", ppi_pair_id)  
        if not os.path.exists(out_desc_dir):  
            os.makedirs(out_desc_dir)  
          
        # 检查是否已经计算过  
        if os.path.exists(os.path.join(out_desc_dir, 'p1_desc_straight.npy')):  
            print("跳过 {} - 描述符已存在".format(ppi_pair_id))  
            continue  
          
        try:  
            tic = time.time()  
            pid = "p1"  
              
            # 加载预计算特征  
            p1_rho_wrt_center = np.load(os.path.join(precomputed_dir, pid + "_rho_wrt_center.npy"))  
            p1_theta_wrt_center = np.load(os.path.join(precomputed_dir, pid + "_theta_wrt_center.npy"))  
            p1_input_feat = np.load(os.path.join(precomputed_dir, pid + "_input_feat.npy"))  
            p1_input_feat = mask_input_feat(p1_input_feat, params["feat_mask"])  
            p1_mask = np.load(os.path.join(precomputed_dir, pid + "_mask.npy"))  
            idx1 = np.array(range(len(p1_rho_wrt_center)))  
              
            print("预计算特征数据加载成功。")  
              
            tic = time.time()  
            print("计算表面描述符 (指纹)...")  
              
            # 计算描述符  
            desc1_str = compute_val_test_desc(  
                learning_obj,  
                idx1,  
                p1_rho_wrt_center,  
                p1_theta_wrt_center,  
                p1_input_feat,  
                p1_mask,  
                batch_size=1000,  
                flip=False,  
            )  
            desc1_flip = compute_val_test_desc(  
                learning_obj,  
                idx1,  
                p1_rho_wrt_center,  
                p1_theta_wrt_center,  
                p1_input_feat,  
                p1_mask,  
                batch_size=1000,  
                flip=True,  
            )  
              
            computation_time = time.time() - tic  
            print("描述符计算完成，耗时: {:.2f}s".format(computation_time))  
              
            # 保存描述符  
            np.save(os.path.join(out_desc_dir, "p1_desc_straight.npy"), desc1_str)  
            np.save(os.path.join(out_desc_dir, "p1_desc_flipped.npy"), desc1_flip)  
              
            print("表面指纹已保存到: {}".format(out_desc_dir))  
            print("Step 4完成! 生成了 {} 个指纹，每个指纹维度为 {}。".format(len(desc1_str), desc1_str.shape[1]))  
              
        except Exception as e:  
            print("处理复合物 {} 时出错: {}".format(ppi_pair_id, str(e)))  
            continue  
      
    print("批量处理完成!")  
  
if __name__ == "__main__":  
    main()