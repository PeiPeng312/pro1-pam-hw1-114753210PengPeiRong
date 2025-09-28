# pro1-pam-hw1-114753210PengPeiRong
NCCU Bioinformatics hw1

可以直接參考與GPT的對話一、二，以下為詳細過程。
對話一：https://chatgpt.com/share/68d930e4-5638-800b-a197-3986847b2463
對話二：https://chatgpt.com/share/68d9333a-4388-800b-90b8-0c2fd00842bb

1.首先請GPT幫我翻譯GITHUB上的README，想先確認題目意涵。（GPT對話內容一連結：https://chatgpt.com/share/68d930e4-5638-800b-a197-3986847b2463）
指令：我有修一堂課叫做"生物資訊概論與實務"，這是作業一，要寫程式，上傳至github及gradescope，我不太知道要怎麼做，請先幫我將老師公布作業的github上的readme翻譯成繁體中文，然後告訴我要怎麼做。(然後貼上README的內容)

GPT翻譯過後告訴我流程為：
(1)用 numpy 或 pandas 讀入 mut.txt。
(2)將突變矩陣數值除以 10,000，確保在 [0,1]。
(3)計算 PAM1 → PAMx（透過矩陣相乘或對角化）。
(4)對每個元素計算 log odds（通常公式是 log(p_ij / (fi * fj))，以 log base 2，再取整數）。
(5)生成一個對稱方陣。
(6)輸出成 pamx.txt，格式與範例一致。

2.接下來我將已有的四份檔案（frequent.png、mut.txt、pam250.txt、hw1_ref.py）丟給GPT，詢問該如何做。
指令：他有提供這四個檔案，我提供給你請你告訴我該怎麼做。（附上三份檔案mut.txt、pam250.txt、hw1_ref.py，以及圖片frequent.png，並附上一張圖片表示他們四個檔案在同一份資料夾內）

GPT告訴我我的流程為：
作業流程（你要做的事）
步驟 1：讀入檔案
步驟 2：計算 PAMx
步驟 3：轉換成對數勝算矩陣 (Log Odds)
步驟 4：輸出檔案
步驟 5：檢查與上傳

3.首先根據方才GPT在對話一給我的流程及程式碼，對照同學提供的步驟，我詢問GPT是否有出入。（GPT對話內容二連結：https://chatgpt.com/share/68d9333a-4388-800b-90b8-0c2fd00842bb）
指令：圖一是同學的提示、步驟。請不要記憶。 以下是你叫我做的內容，請對比一下是否正確。（然後貼上對話一的步驟及程式碼）

GPT發現不同之處如下：
📌 不同的地方（你要修正的）
一、對數底數不同
你的版本寫 log₂，但同學的提示（老師給的公式）是
．所以要用10 x log10，不是log₂。
．而且分母只有fi，不是fi x fi。
二、歸一化 PAM1 的方式
．同學的步驟是「mut/10000 得到機率 M1」，再直接做矩陣乘方。
．你的版本有多一步「行內正規化 (row sum = 1)」。
．如果 mut.txt 本來就是設計成「次數/10000 = 機率」，那可能不需要再正規化。
．建議你照老師的步驟：直接除以 10000 得到機率。

GPT提供的修正要點如下：
📌 修正版程式要點
．10 x log10，不是log₂。
．分母用單一頻率fi​，不是乘積。
．M1 = mut/10000，不需要再 row-normalize。

4.因為GPT提供的Log-Odds Score公式不對，所以回去對話1貼了簡報上的公式給它（sij​=10×log10​Rij​,其中 Rij​=​Mij​​/fi）（ij都是下標）。
　而對話二我也再次貼上簡報上的Log-Odds Score公式，它告訴我下一步驟是「把 frequent.png 的 各胺基酸頻率數值抄成一個字典」。

5.指令：請一步驟一步驟地告訴我程式碼。
步驟 1：建立 function
步驟 2：將 frequent.png 的頻率寫成常數
步驟 3 & 4：讀取 mut.txt，並除以 10000 得到 M1
步驟 5：M1 自乘 x 次，得到 Mx
步驟 6：計算 Log-Odds Score
步驟 7 & 8：輸出 pamx.txt

6.由於看不懂frequent.png胺基酸頻率是如何定義的因此詢問GPT緣由，詢問後並對應維基百科的"標準蛋白胺基酸列表"去檢查是否正確。

7.接著請一步驟一步驟地問GPT該怎麼寫Python並檢查是否有誤，每步驟都會叫它整合當前所有步驟。

8.在整合完上述步驟1~8程式碼後，GPT建議我在程式碼中加上「if __name__ == "__main__":」，我請它幫我加入，它就幫我輸出完整的python程式碼。
指令：好的請幫我加上if __name__ == "__main__"

9.在執行程式的時候發現Error，於是我將結果及Error內容丟給GPT詢問問題在哪。
指令：結果出現了 === 次數矩陣 (mut) === [] === 機率矩陣 (M1) === [] 以及error：LinAlgError: Last 2 dimensions of the array must be square

GPT表示：
📌 問題出在 mut.txt 沒有正確讀進來
可能原因有三種：
(1)檔案路徑錯誤：mut.txt 不在你執行程式的工作目錄。
(2)檔案分隔符號不是 \t (tab)：可能是空白格或逗號。
(3)檔案格式有多餘標頭：第一列不是你想要的標題。

📌 解法2：查看mut.txt前五行
with open("mut.txt", "r") as f:
    for i in range(5):
        print(f.readline())
        
10.根據建議我檢查mut.txt檔案前幾行，對照GPT剛剛的程式碼發現有誤，原先程式碼以為mut.txt是用tab（也就是\t）當分隔符號，於是我將mut.txt及查看前5行的程式碼的結果貼給GPT，發現確實是讀取mut.txt時出錯。更正前後的讀取mut.txt程式碼如下：
更正前：
pd.read_csv(input_path, sep="\t", index_col=0)
更正後：
df = pd.read_csv(
    input_path,
    sep=r"\s+",        # 用正則表達式，自動處理多個空白
    comment="#",        # 跳過註解行
    index_col=0         # 第一欄 (A, R, N...) 當 row label
)

11.最終請GPT將更正後的程式碼結合前面的程式碼。
指令：好的請幫我將修正版結合進完整的程式碼再給我一次程式碼。

12.執行完程式碼後確認成功輸出pamx.txt後，請GPT給我程式碼讓我比對答案是否跟pam250.txt相同。
指令：我要怎麼知道答案是否正確。

13.最後在程式碼最前方輸入學號、姓名，以及將程式碼的輸出檔案名改成pamx.txt。
