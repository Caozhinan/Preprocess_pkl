import numpy as np    
import os    
import sys    
import time    
import importlib    
    
# 导入必要的模块    
from default_config.masif_opts import masif_opts    
from masif_modules.MaSIF_ppi_search import MaSIF_ppi_search    
from masif_modules.train_ppi_search import compute_val_test_desc    
    
def mask_input_feat(input_feat, mask):    
    """    
    根据特征掩码过滤输入特征。    
    """    
    mymask = np.where(np.array(mask) == 0.0)[0]    
    return np.delete(input_feat, mymask, axis=2)    
    
def generate_surface_fingerprints(precomputed_features_dir, output_dir, ppi_pair_id, custom_params_file=None):    
    """    
    Step 4: 基于神经网络的指纹生成    
        
    参数:    
    - precomputed_features_dir: Step 3生成的预计算特征的目录    
    - output_dir: 最终指纹的输出目录    
    - ppi_pair_id: 蛋白质-配对ID (例如 'complex_A')    
    - custom_params_file: 包含自定义参数的Python文件路径 (例如 'nn_models.sc05.all_feat.custom_params')    
        
    返回:    
    - 包含指纹数据的字典    
    """    
        
    print(f"开始Step 4: 基于神经网络的指纹生成...")    
    # print(f"从目录加载预计算特征: {precomputed_features_dir}")    
        
    # 1. 设置参数    
    params = masif_opts["ppi_search"].copy()  # 创建副本避免修改原始配置  
    if custom_params_file:      
        try:  
            # 检查文件是否存在  
            if os.path.exists(custom_params_file):  
                # 将文件路径转换为模块路径格式  
                module_path = custom_params_file.replace('.py', '').replace('/', '.')  
                # 如果是绝对路径，需要添加到sys.path  
                import sys  
                file_dir = os.path.dirname(custom_params_file)  
                if file_dir not in sys.path:  
                    sys.path.insert(0, file_dir)  

                # 只使用文件名作为模块名  
                module_name = os.path.basename(custom_params_file).replace('.py', '')  
                custom_params_module = importlib.import_module(module_name)  
                custom_params = custom_params_module.custom_params  

                for key in custom_params:  
                    # print("设置 {} 为 {}".format(key, custom_params[key]))  
                    params[key] = custom_params[key]  
            else:  
                print("自定义参数文件不存在: {}".format(custom_params_file))  
                print("将使用默认参数。")  
        except Exception as e:  
            print("警告: 无法加载自定义参数文件 {}: {}".format(custom_params_file, str(e)))  
            # print("将使用默认参数。")
    # 2. 加载预训练的神经网络模型    
    # print("加载预训练的神经网络模型...")    
    learning_obj = MaSIF_ppi_search(    
        params["max_distance"],    
        n_thetas=16,    
        n_rhos=5,    
        n_rotations=16,    
        idx_gpu="/gpu:0",    
        feat_mask=params["feat_mask"],    
    )    
      
# 恢复模型权重 - 改进的路径处理  
# 首先检查 HDF5 权重文件（绝对路径）  
    hdf5_weight_dir = "/es01/paratera/sce0413/czn/preprocess_pkl/masif/data/masif_ppi_search/nn_models/sc05/all_feat/weight"  
    hdf5_weight_path = os.path.join(hdf5_weight_dir, "weights_12A_0129.hdf5")  
    # model_load_start = time.time()
    if os.path.exists(hdf5_weight_path + ".data-00000-of-00001"):  
        model_path = hdf5_weight_path  
        # print(f"使用HDF5权重文件: {hdf5_weight_path}")  
    else:  
        # 回退到原始的 checkpoint 路径检查  
        model_path = os.path.join(params["model_dir"], "model")  
          
        if not os.path.exists(model_path + ".meta"):  
            # 尝试绝对路径  
            abs_model_dir = "/es01/paratera/sce0413/czn/preprocess_pkl/masif/data/masif_ppi_search/nn_models/sc05/all_feat/model_data/"  
            abs_model_path = os.path.join(abs_model_dir, "model")  
              
            if os.path.exists(abs_model_path + ".meta"):  
                model_path = abs_model_path  
                # print(f"使用绝对路径加载模型: {abs_model_dir}")  
            else:  
                print("错误: 未找到可用的模型文件")  
                return None  
        else:  
            pass
        try:  
            learning_obj.saver.restore(learning_obj.session, model_path)    
            # print("神经网络模型加载成功。")  
        except Exception as e:  
            print(f"模型加载失败: {e}")  
            return None  
    # model_load_time = time.time() - model_load_start
    # print(f"[fingerprint_gen.py] 模型权重加载耗时: {model_load_time:.2f} 秒")
    # 3. 加载预计算的特征数据    
    in_dir = precomputed_features_dir    
        
    try:      
        # 从压缩文件加载所有特征  
        npz_file = np.load(os.path.join(in_dir, "all_features.npz"))  
          
        p1_rho_wrt_center = npz_file['p1_rho_wrt_center']  
        p1_theta_wrt_center = npz_file['p1_theta_wrt_center']  
        p1_input_feat = npz_file['p1_input_feat']  
        p1_mask = npz_file['p1_mask']  
          
        npz_file.close()  # 释放文件句柄  
              
        # 过滤输入特征 (如果feat_mask有定义)      
        p1_input_feat = mask_input_feat(p1_input_feat, params["feat_mask"])      
              
        idx1 = np.array(range(len(p1_rho_wrt_center)))      
    except Exception as e:      
        print(f"错误: 无法加载预计算特征数据: {e}")      
        print(f"请确保文件位于 {in_dir} 且命名正确。")      
        return None 
    
    # 4. 计算表面描述符 (指纹)    
    # print("计算表面描述符 (指纹)...")    
    # tic = time.time()    
        
    # 计算标准方向的描述符    
    desc1_str = compute_val_test_desc(    
        learning_obj,    
        idx1,    
        p1_rho_wrt_center,    
        p1_theta_wrt_center,    
        p1_input_feat,    
        p1_mask,    
        batch_size=24, # 可以根据GPU内存调整    
        flip=False,    
    )    
        
    # 计算翻转方向的描述符 (用于互补性分析)    
    desc1_flip = compute_val_test_desc(    
        learning_obj,    
        idx1,    
        p1_rho_wrt_center,    
        p1_theta_wrt_center,    
        p1_input_feat,    
        p1_mask,    
        batch_size=24, # 可以根据GPU内存调整    
        flip=True,    
    )    
    # print(f"描述符计算完成，耗时: {time.time() - tic:.2f}s")    
    
    # 5. 保存描述符    
    out_desc_dir = os.path.join(output_dir, "descriptors", ppi_pair_id)    
    os.makedirs(out_desc_dir, exist_ok=True)    
        
    np.save(os.path.join(out_desc_dir, "p1_desc_straight.npy"), desc1_str)    
    np.save(os.path.join(out_desc_dir, "p1_desc_flipped.npy"), desc1_flip)    
    # print(f"表面指纹已保存到: {out_desc_dir}")    
    
    # 6. 构建返回的指纹字典    
    fingerprints = {    
        'desc_straight': desc1_str,    
        'desc_flipped': desc1_flip,    
        'output_dir': out_desc_dir    
    }    
        
    print(f"Step 4完成! 生成了 {desc1_str.shape[0]} 个指纹，每个指纹维度为 {desc1_str.shape[1]}。")    
        
    return fingerprints    
    
# 使用示例    
if __name__ == "__main__":
    import argparse
    start_time = time.time()
    parser = argparse.ArgumentParser(description='生成MaSIF表面指纹')
    parser.add_argument('--precomputed_dir', required=True, help='预计算特征目录')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    parser.add_argument('--ppi_pair_id', required=True, help='蛋白质对ID')
    parser.add_argument('--custom_params_file', default="masif.source.masif_ppi_search.nn_models.sc05.all_feat.custom_params",
                       help='自定义参数文件')

    args = parser.parse_args()

    fingerprints = generate_surface_fingerprints(
        precomputed_features_dir=args.precomputed_dir,
        output_dir=args.output_dir,
        ppi_pair_id=args.ppi_pair_id,
        custom_params_file=args.custom_params_file
    )

    if fingerprints is not None:
        # print("\n=== Step 4 输出指纹总结 ===")
        print(f"标准方向指纹形状: {fingerprints['desc_straight'].shape}")
        # print(f"翻转方向指纹形状: {fingerprints['desc_flipped'].shape}")
        # print(f"指纹文件保存在: {fingerprints['output_dir']}")
    else:
        print("指纹生成失败!")
        sys.exit(1)

    elapsed = time.time() - start_time
    print(f"[fingerprint_gen.py] Finished in {elapsed:.2f} seconds")