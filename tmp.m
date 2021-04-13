load('../loss_save_all')

lambda_list = [0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2];
shift_list = [0, 5, 10, 15];

imagesc(x / max(max(x)))
colorbar()
