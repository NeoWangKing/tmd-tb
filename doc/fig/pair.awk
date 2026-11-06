# 把两个文件的数据行并排拼起来（跳过 # 注释行与空行），供 gnuplot 算列间差
# 用法: awk -f doc/fig/pair.awk 文件1 文件2
FNR == NR {
    if ($1 !~ /^#/ && NF > 2) a[++n] = $0
    next
}
$1 !~ /^#/ && NF > 2 { print a[++m], $0 }
