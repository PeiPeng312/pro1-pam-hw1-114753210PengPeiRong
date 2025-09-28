# hw1.py
# Name: 彭珮蓉
# Student ID: 114753210

import numpy as np
import pandas as pd

def generate_pam(x, input_path, output_path, freqs=None):
    # Step1：胺基酸頻率表 (frequent.png) (對照了 維基百科-標準蛋白胺基酸列表)
    if freqs is None:
        freqs = {
            "G": 0.089, "A": 0.087, "L": 0.085, "K": 0.081,
            "S": 0.070, "V": 0.065, "T": 0.058, "P": 0.051,
            "E": 0.050, "D": 0.047, "R": 0.041, "N": 0.040,
            "F": 0.040, "Q": 0.038, "I": 0.037, "H": 0.034,
            "C": 0.033, "Y": 0.030, "M": 0.015, "W": 0.010
        }

    # Step2：讀取突變次數矩陣 (mut.txt)
    df = pd.read_csv(
        input_path,
        sep=r"\s+",        # 空白分隔
        comment="#",       # 跳過註解行
        index_col=0        # 第一欄 (A,R,N...) 當 row label
    )
    mut_matrix = df.values
    amino_acids = df.index.tolist()

    print("=== 次數矩陣 (mut) ===")
    print(mut_matrix[:5])

    # Step3：除以 10000 → 機率矩陣 M1
    M1 = mut_matrix / 10000.0

    print("\n=== 機率矩陣 (M1) ===")
    print(M1[:5])

    # Step4：M1 自乘 x 次 → Mx
    Mx = np.linalg.matrix_power(M1, x)

    print(f"\n=== M{x} 矩陣 ===")
    print(Mx[:5])

    # Step5：計算 Log-Odds Score → PAMx
    PAMx = np.zeros_like(Mx, dtype=int)
    for i, aa_i in enumerate(amino_acids):
        for j, aa_j in enumerate(amino_acids):
            Rij = Mx[i, j] / freqs[aa_i]   # Mij / fi
            score = 10 * np.log10(Rij)     # 10 * log10
            PAMx[i, j] = int(round(score)) # 四捨五入

    print(f"\n=== PAM{x} 矩陣 (整數分數) ===")
    print(PAMx[:5])

    # Step6：輸出 pamx.txt 檔案
    with open(output_path, "w") as f:
        f.write("\t" + "\t".join(amino_acids) + "\n")
        for i, aa in enumerate(amino_acids):
            row = "\t".join(map(str, PAMx[i]))
            f.write(aa + "\t" + row + "\n")

    print(f"\n 已輸出結果到 {output_path}")

    return PAMx, amino_acids


# === 主程式入口 ===
if __name__ == "__main__":
    # 預設跑 PAM250
    PAM250, amino_acids = generate_pam(
        250,               # x = PAM250
        "mut.txt",         # 輸入檔 (與 hw1.py 同資料夾)
        "pamx.txt"   # 輸出檔
    )