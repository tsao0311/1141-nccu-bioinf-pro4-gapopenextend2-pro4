import sys
from Bio import SeqIO
import numpy as np

# 使用一個極大負數作為不可能路徑的得分
IMPOSSIBLE_SCORE = -1e20 

# --- 步驟 1: 讀取替換得分矩陣 ---
def load_score_matrix(score_path):
    """
    從檔案載入替換得分矩陣 (e.g., PAM250)。
    返回一個字典，鍵為(char1, char2)，值為得分。
    """
    matrix = {}
    header = []
    
    # 設置默認值以處理未知字符
    default_score = -100 # 假設的極大懲罰
    
    try:
        with open(score_path, 'r') as f:
            lines = f.readlines()
            
            # 找到標頭行
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if not header:
                    # 設置標頭，排除第一個可能是空或註釋符號的元素
                    header = [p for p in parts if len(p) == 1]
                    continue
                
                # 處理數據行
                row_char = parts[0]
                # 確保只處理包含得分的行，且行標識是有效的氨基酸
                if row_char in header: 
                    # 從第二個元素開始是分數
                    scores = [int(s) for s in parts[1:]] 
                    
                    for j, col_char in enumerate(header):
                        if j < len(scores):
                             # 儲存得分，考慮矩陣對稱性
                            matrix[(row_char, col_char)] = float(scores[j])
                            matrix[(col_char, row_char)] = float(scores[j])
                            
    except FileNotFoundError:
        print(f"錯誤：找不到得分檔案 {score_path}")
        sys.exit(1)
    except Exception as e:
        print(f"讀取得分檔案時發生錯誤: {e}")
        sys.exit(1)
        
    return matrix, header

# --- 輔助函數: 寫入 FASTA 檔案 ---
def write_fasta(output_path, seq_id_1, aln_seq_1, seq_id_2, aln_seq_2, score):
    """將對齊結果寫入 FASTA 格式的檔案。"""
    try:
        # 將分數格式化為一位小數，以匹配您的測試輸出
        score_str = f"{score:.1f}" 
        with open(output_path, 'w') as f:
            f.write(f">{seq_id_1} Score: {score_str}\n") 
            f.write(f"{aln_seq_1}\n")
            f.write(f">{seq_id_2} Score: {score_str}\n")
            f.write(f"{aln_seq_2}\n")
    except Exception as e:
        print(f"寫入輸出檔案時發生錯誤: {e}")
        
# --- 步驟 3 & 4: 實現 Gotoh 動態規劃與回溯 (M/Ix/Iy 邏輯最終修正) ---
def gotoh_alignment(seq1, seq2, score_matrix, gap_open, gap_extend, aln_type):
    """
    使用 Gotoh 演算法進行序列對齊 (Needleman-Wunsch/Smith-Waterman)。
    """
    n, m = len(seq1), len(seq2)
    
    # 初始化 DP 矩陣：使用 float64 確保精度
    M = np.zeros((n + 1, m + 1), dtype=np.float64)
    Ix = np.full((n + 1, m + 1), IMPOSSIBLE_SCORE, dtype=np.float64) 
    Iy = np.full((n + 1, m + 1), IMPOSSIBLE_SCORE, dtype=np.float64) 
    
    # 回溯矩陣 (TB = Traceback)
    TB_M = np.zeros((n + 1, m + 1), dtype=int)
    TB_Ix = np.zeros((n + 1, m + 1), dtype=int)
    TB_Iy = np.zeros((n + 1, m + 1), dtype=int)

    # --- 邊界初始化 (全域對齊) ---
    for i in range(1, n + 1):
        score = gap_open + i * gap_extend 
        M[i, 0] = score
        Ix[i, 0] = score # 必須由 Ix 處理 Seq1 的間隙
        # Iy[i, 0] 保持 IMPOSSIBLE_SCORE
        TB_M[i, 0] = 1 
        TB_Ix[i, 0] = 1 
        
    for j in range(1, m + 1):
        score = gap_open + j * gap_extend
        M[0, j] = score
        Iy[0, j] = score # 必須由 Iy 處理 Seq2 的間隙
        # Ix[0, j] 保持 IMPOSSIBLE_SCORE
        TB_M[0, j] = 2 
        TB_Iy[0, j] = 1 

    # 局部對齊 (Smith-Waterman) 的特殊邊界初始化：
    if aln_type == 'local':
        # 在局部對齊中，所有負值都會被 0 取代，所以只需確保 M[0,0] 為 0
        M[:, :] = 0.0 
        Ix[:, :] = 0.0
        Iy[:, :] = 0.0
        # 重設回溯矩陣
        TB_M[:, :] = 0 
        TB_Ix[:, :] = 0
        TB_Iy[:, :] = 0
        
    # 動態規劃填充
    max_score = 0.0
    max_i, max_j = n, m
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            
            # 1. 獲取替換得分 (S)。如果找不到，假設其懲罰極高
            s = score_matrix.get((seq1[i-1], seq2[j-1]), IMPOSSIBLE_SCORE) 
            
            # 2. M 矩陣的遞歸關係 (Match/Mismatch)
            # M[i, j] = S(A_i, B_j) + max(M[i-1, j-1], Ix[i-1, j-1], Iy[i-1, j-1])
            prev_M = M[i-1, j-1]
            prev_Ix = Ix[i-1, j-1]
            prev_Iy = Iy[i-1, j-1]

            # 計算三個來源路徑的分數 (這裡的 s 應用於所有來源)
            score_M = prev_M + s
            score_Ix = prev_Ix + s
            score_Iy = prev_Iy + s
            
            # 選擇 M 矩陣的最佳前導路徑
            M[i, j] = max(score_M, score_Ix, score_Iy)
            
            # 記錄 M 的回溯路徑
            if M[i, j] == score_M:
                TB_M[i, j] = 0 # 來自 M[i-1, j-1]
            elif M[i, j] == score_Ix:
                TB_M[i, j] = 1 # 來自 Ix[i-1, j-1]
            else:
                TB_M[i, j] = 2 # 來自 Iy[i-1, j-1]
                
            # 3. Ix 矩陣 (Gap in Seq2, Seq1 向下延伸)
            # 開啟間隙：來自 M[i-1, j] + G_open + G_extend
            open_gap_ix = M[i-1, j] + gap_open + gap_extend
            # 延伸間隙：來自 Ix[i-1, j] + G_extend
            extend_gap_ix = Ix[i-1, j] + gap_extend
            Ix[i, j] = max(open_gap_ix, extend_gap_ix)
            
            # 記錄 Ix 的回溯路徑
            if Ix[i, j] == open_gap_ix:
                TB_Ix[i, j] = 0 # 來自 M[i-1, j]
            else:
                TB_Ix[i, j] = 1 # 來自 Ix[i-1, j]
                
            # 4. Iy 矩陣 (Gap in Seq1, Seq2 向右延伸)
            # 開啟間隙：來自 M[i, j-1] + G_open + G_extend
            open_gap_iy = M[i, j-1] + gap_open + gap_extend
            # 延伸間隙：來自 Iy[i, j-1] + G_extend
            extend_gap_iy = Iy[i, j-1] + gap_extend
            Iy[i, j] = max(open_gap_iy, extend_gap_iy)
            
            # 記錄 Iy 的回溯路徑
            if Iy[i, j] == open_gap_iy:
                TB_Iy[i, j] = 0 # 來自 M[i, j-1]
            else:
                TB_Iy[i, j] = 1 # 來自 Iy[i, j-1]

            # 局部對齊 (Smith-Waterman) 的額外條件
            if aln_type == 'local':
                # 所有負值皆重置為 0
                M[i, j] = max(0.0, M[i, j])
                Ix[i, j] = max(0.0, Ix[i, j])
                Iy[i, j] = max(0.0, Iy[i, j])
                
                # 更新最高得分和其位置 (滿足最長對齊的簡化版 tie-breaker)
                current_max = max(M[i, j], Ix[i, j], Iy[i, j])
                if current_max >= max_score: 
                    max_score = current_max
                    max_i, max_j = i, j

    # --- 決定最終得分和起始回溯點 ---
    if aln_type == 'global':
        final_score = max(M[n, m], Ix[n, m], Iy[n, m])
        # 設置回溯的起始矩陣 (優先 M > Ix > Iy，以處理平分情況)
        if M[n, m] == final_score:
            current_matrix = 0 
        elif Ix[n, m] == final_score:
            current_matrix = 1
        else:
            current_matrix = 2
        start_i, start_j = n, m

    elif aln_type == 'local':
        final_score = max_score
        # 設置回溯的起始矩陣 (在 max_i, max_j 處)
        current_max_at_max_pos = max(M[max_i, max_j], Ix[max_i, max_j], Iy[max_i, max_j])
        if current_max_at_max_pos == M[max_i, max_j]:
            current_matrix = 0
        elif current_max_at_max_pos == Ix[max_i, max_j]:
            current_matrix = 1
        else:
            current_matrix = 2
        start_i, start_j = max_i, max_j
    
    # --- 回溯 ---
    aln_seq1 = ""
    aln_seq2 = ""
    i, j = start_i, start_j
    
    while i > 0 or j > 0:
        if aln_type == 'local' and max(M[i, j], Ix[i, j], Iy[i, j]) <= 0:
             # 局部對齊遇到 0 分或負分時停止
             break 
        
        if current_matrix == 0: # M 矩陣 (比對)
            trace = TB_M[i, j]
            aln_seq1 = seq1[i-1] + aln_seq1
            aln_seq2 = seq2[j-1] + aln_seq2
            i -= 1
            j -= 1
            current_matrix = trace # TB_M[i,j] 記錄的是下一個矩陣 (0, 1, or 2)
                
        elif current_matrix == 1: # Ix 矩陣 (Seq2 有間隙，Seq1 向下延伸)
            trace = TB_Ix[i, j]
            aln_seq1 = seq1[i-1] + aln_seq1
            aln_seq2 = '-' + aln_seq2
            i -= 1
            current_matrix = 0 if trace == 0 else 1
                
        elif current_matrix == 2: # Iy 矩陣 (Seq1 有間隙，Seq2 向右延伸)
            trace = TB_Iy[i, j]
            aln_seq1 = '-' + aln_seq1
            aln_seq2 = seq2[j-1] + aln_seq2
            j -= 1
            current_matrix = 0 if trace == 0 else 2

    return final_score, aln_seq1, aln_seq2


# --- 步驟 5: 主函數 (保持不變) ---
def alignment(input_path, score_path, output_path, aln, gap_open, gap_extend):
    """
    執行成對序列比對的主函數。
    """
    
    # 1. 載入得分矩陣
    score_matrix, _ = load_score_matrix(score_path)

    # 2. 解析 FASTA 檔案並讀取序列
    sequences = []
    ids = []
    try:
        for record in SeqIO.parse(input_path, "fasta"):
            sequences.append(str(record.seq).upper()) # 確保序列為大寫
            ids.append(record.id)
            if len(sequences) == 2:
                break
    except Exception as e:
        print(f"解析 FASTA 檔案時發生錯誤: {e}")
        return

    if len(sequences) < 2:
        print("錯誤：FASTA 檔案中需要至少兩個序列進行成對比對。")
        return
        
    seq1, seq2 = sequences[0], sequences[1]
    id1, id2 = ids[0], ids[1]
    
    # 3. 執行 Gotoh 對齊
    final_score, aln_seq1, aln_seq2 = gotoh_alignment(
        seq1, seq2, score_matrix, gap_open, gap_extend, aln
    )

    # 4. 寫入結果
    write_fasta(output_path, id1, aln_seq1, id2, aln_seq2, final_score)
    
    # 輸出結果格式與您提供的測試輸出一致
    print(f"對齊完成。得分: {final_score:.1f}，結果已寫入 {output_path}")

# --- 執行入口 (保持不變) ---
if __name__ == "__main__":
    
    # 請根據您的腳本調用方式，在此處添加命令行解析邏輯或直接調用 alignment 函數
    print("請使用您的命令行參數調用 alignment 函數以進行測試。")
