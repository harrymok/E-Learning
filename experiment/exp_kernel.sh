exp=${1:-1}
case=${2:-1}
id=${3:-1}

case $case in
1)
  tf=1
  prop=1
  hetero=0
  ;;
2)
  tf=0
  prop=1
  hetero=0
  ;;
3)
  tf=1
  prop=1
  hetero=1
  ;;
4)
  tf=0
  prop=1
  hetero=1
  ;;
5)
  tf=1
  prop=0
  hetero=0
  ;;
6)
  tf=0
  prop=0
  hetero=0
  ;;
7)
  tf=1
  prop=0
  hetero=1
  ;;
8)
  tf=0
  prop=0
  hetero=1
  ;;
*)
  echo "incorrect case number: $case"
  exit 1
esac

n=400
p=10
K=3

case $exp in
1)
  ;;
# n in 100 200 800 1600
2)
  n=100
  ;;
3)
  n=200
  ;;
4)
  n=800
  ;;
5)
  n=1600
  ;;
# K = 2
6)
  K=2
  ;;
7)
  n=100
  K=2
  ;;
8)
  n=200
  K=2
  ;;
9)
  n=800
  K=2
  ;;
10)
  n=1600
  K=2
  ;;
# p in 50 100
11)
  p=50
  ;;
12)
  n=100
  p=50
  ;;
13)
  n=200
  p=50
  ;;
14)
  n=800
  p=50
  ;;
15)
  n=1600
  p=50
  ;;
16)
  p=100
  ;;
17)
  n=100
  p=100
  ;;
18)
  n=200
  p=100
  ;;
19)
  n=800
  p=100
  ;;
20)
  n=1600
  p=100
  ;;
# K in 5 7
21)
  K=5
  ;;
22)
  n=100
  K=5
  ;;
23)
  n=200
  K=5
  ;;
24)
  n=800
  K=5
  ;;
25)
  n=1600
  K=5
  ;;
26)
  K=7
  ;;
27)
  n=100
  K=7
  ;;
28)
  n=200
  K=7
  ;;
29)
  n=800
  K=7
  ;;
30)
  n=1600
  K=7
  ;;
# p = 50, K = 2
31)
  p=50
  K=2
  ;;
32)
  n=100
  p=50
  K=2
  ;;
33)
  n=200
  p=50
  K=2
  ;;
34)
  n=800
  p=50
  K=2
  ;;
35)
  n=1600
  K=2
  p=50
  ;;
*)
  echo "incorrect exp number: $exp"
  exit 1
esac

dir0="RData"
dir01="${dir0}/exp_${exp}"
dir02="${dir01}/case_${case}"
save="${dir02}/exp_${id}.RData"

dir1="Rout"
out="${dir1}/out_${exp}_${case}_${id}.Rout"

[ ! -d $dir0 ] && mkdir $dir0
[ ! -d $dir01 ] && mkdir $dir01
[ ! -d $dir02 ] && mkdir $dir02
[ ! -d $dir1 ] && mkdir $dir1

setup="
  dir='..'
  file_save='$save' 

  n=$n
  p=$p
  K=$K
  treatment_free_effect=$tf
  propensity_score=$prop
  heteroscedasticity=$hetero
"

if [ -f $save ]
then
  echo "override $save"
  rm $save
fi
R CMD BATCH --no-save --no-restore "--args $setup" exp_kernel.R $out
