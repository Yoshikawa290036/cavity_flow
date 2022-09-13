# cavity_flow
2次元キャビティー流れのシミュレーション

Linuxで，OpenMP，g++，gnuplot，Image-Magick，texlive，ffmpeg がインストールされている環境下で行う．
また，Image-Magickは，eps → png に convert できるよう，policy.xmlを編集しておく．

1. "data"というディレクトリを作成する．

$ mkdir data


2. 数値計算する cavity.cpp　をコンパイルする．

$ g++ cavity.cpp -fopenmp -O3


3. コンパイル後，生成されたフィアル (a.out)を実行する．すると，"data"ディレクトリにデータファイルが作成される．

$ ./a.out


4. 3.で生成された流れ場を可視化する．call_mk_img.shを実行すると，mk_img.sh が呼び出され，"image"ディレクトリ内に，可視化された流れ場の画像が生成される．

$ ./call_mk_img.sh


5. imageディレクトリ内の画像から動画を作成する．mk_mp4.sh を実行すると cavity_flow.mp4 のような動画が作成される．

$ ./mk_mp4.sh
