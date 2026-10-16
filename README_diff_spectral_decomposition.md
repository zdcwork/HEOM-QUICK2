# 不同电极谱分解方案说明

## 1. 修改目标

本次修改在尽量少改动原 HEOM 程序结构的前提下，新增 diff 模式，使不同环境电子库/不同电极 `ialf` 可以使用不同谱分解方案。原始 `main/evaluatepara.f90` 保持不变；diff 模式由新增的 `main/evaluatepara_diff.f90` 负责。

## 2. 涉及文件

- `include/common`
  - 新增 `$diffdecomp` 相关全局变量和逐电极配置数组。
- `main/main.f90`
  - 预读 `$diffdecomp` 开关。
  - 普通模式继续调用 `evaluatepara(init)`。
  - diff 模式调用 `evaluatepara_diff(init)`，并跳过原全局 scheme 互斥检查。
- `main/evaluatepara_diff.f90`
  - 从原 `evaluatepara.f90` 复制而来，子程序名改为 `evaluatepara_diff`。
  - 新增逐电极 `cb/cd/cgama` 分支。
- `main/read_diffdecomp_config.f90`
  - 读取 `$diffdecomp` 和后续 `nalf` 行 scheme 配置。
- `main/diff_get_psf_params.f90`
  - 按 `itype_psf` 返回 PSF 预设 pole 数、`rpsf`、`tpsf`。
- `tools/read_diff_check/read_diff_check.f90`
  - 独立 parser 测试程序源码，不在 `main/*.f90` 中，避免产生第二个主程序。
- `README_diff_spectral_decomposition.md`
  - 本说明文档。

## 3. 新增变量

新增变量位于 `include/common`：

```fortran
logical :: ldiff_spe_dec, ldiff_spe_dec_sp
integer :: ncount_alf(maxalf)
integer :: ncount_scheme(maxalf)
integer :: ncount_poles(maxalf)
integer :: ncount_type(maxalf)
integer :: nmats_alf(maxalf)
integer :: ncor_alf(maxalf)
logical :: lpsdfff_alf(maxalf)
```

说明：

- `ldiff_spe_dec/ldiff_spe_dec_sp`：是否启用逐电极谱分解。
- `ncount_alf(ialf)`：配置中指定的电极编号。
- `ncount_scheme(ialf)`：该电极 scheme 编码。
- `ncount_poles(ialf)`：输入 pole 数。
- `ncount_type(ialf)`：Fano/PSF 的 type 参数。
- `nmats_alf(ialf)`：该电极真实谱展开项数。
- `ncor_alf(ialf)`：该电极真实参与 HEOM 的 correlation 项数。
- `lpsdfff_alf(ialf)`：该电极是否为 Fano/PSD-FFF。

为保持 common-block 风格和 Makefile 编译顺序，数组采用 `maxalf` 固定维度，而不是 Fortran module allocatable 数组。

## 4. input 文件格式

普通模式可不写 `$diffdecomp`，或写：

```fortran
$diffdecomp ldiff_spe_dec = .false. ldiff_spe_dec_sp = .false. $end
```

逐电极 diff 模式示例：

```fortran
$diffdecomp ldiff_spe_dec = .true. ldiff_spe_dec_sp = .false. $end
1   Pade   20
2   fano   3  3
```

每行格式：

```text
ialf   scheme_name   n_poles_or_nmats   optional_type
```

`optional_type` 对 Fano 和 PSF 必须为正整数。

## 5. scheme 字典

```text
1  pade    -> Pade/PSD
2  psdjob  -> Pade/PSD
3  prony   -> PSF/预设 Prony
4  psfjob  -> PSF/预设 Prony
5  psdfff  -> Fano/PSD-FFF
6  fano    -> Fano/PSD-FFF
7  read    -> 读取外部 Prony 分解结果
8  mtsjob  -> Matsubara
9  lprony  -> 读取外部 Prony 分解结果
```

scheme 名称大小写不敏感，例如 `Pade`、`PSDJob`、`Fano`、`MTSJOB` 都可识别。

## 6. main.f90 调用逻辑

`main.f90` 预读 `$diffdecomp` 的两个开关：

```fortran
if (ldiff_spe_dec .or. ldiff_spe_dec_sp) then
    call evaluatepara_diff(init)
else
    call evaluatepara(init)
end if
```

普通模式仍执行原来的全局 scheme 互斥检查。diff 模式跳过该全局互斥检查，由 `read_diffdecomp_config` 和 `evaluatepara_diff` 按电极检查。

## 7. evaluatepara_diff 核心实现

`evaluatepara_diff` 读完 input 基本参数并得到 `nalf` 后，调用：

```fortran
call read_diffdecomp_config(5, nalf)
```

随后按电极设置：

- Pade/Fano/Matsubara/read：`nmats_alf(ialf)=ncount_poles(ialf)`。
- PSF：由 `ncount_type(ialf)` 查询预设 `npsf/rpsf/tpsf`。若 input pole 数与预设不一致，输出 warning，并以预设 pole 数为准。
- 全局 `nmats=maxval(nmats_alf)`。
- 全局 `ncor=maxval(ncor_alf)`。

`cb/cd/cgama` 填充时，在 `do nj=1,nalf` 内按 `ncount_scheme(nj)` 分流：

- Pade：复用原 PSD 公式。
- Fano：复用原 PSD-FFF 公式。
- PSF：使用该电极自己的 `cpcoef_psf_alf(:,nj)`、`cppole_psf_alf(:,:,nj)`、`tpsf_alf(nj)`。
- Matsubara：显式使用原 Matsubara 公式。
- read/lprony：只对该电极从 `PRONY_CB/CD/GAMA-ISGN-*-IALF-*-IORBS-*.data` 读取。

## 8. padding 策略

由于旧程序数组仍使用统一全局 `ncor`，短展开电极需要 padding。对无真实贡献项：

```fortran
cb(...)    = czero
cd(...)    = czero
cgama(...) = dcmplx(-1.d30, 0.d0)
```

这样可以避免 `cgama=0` 导致后续除零或 NaN。`dimpf` 只在真实项且 `abs(cgama)>dpico` 时计算。

Fano 额外项统一放在全局 `ndrude + nmats + 1` 之后，以兼容原有 `mfff1d/mmfff/ctheta` 抽屉标记逻辑。某个 Fano 电极若 `nmats_alf < nmats`，中间空位同样用 padding。

## 9. 编译和测试

原工程编译方式不变：

```bash
make
```

当前环境限制：

- 当前 PowerShell 环境没有 `make`。
- 当前环境没有 `ifort`。
- 使用 `gfortran -fsyntax-only` 对新增小文件和 `evaluatepara_diff.f90` 做了语法检查，均通过。
- `main.f90` 的 gfortran 检查被原有 `module/sj_timer.f90` 中 Intel 风格格式串阻挡，不是本次修改引入的问题。

建议测试：

1. 普通模式：

```fortran
$diffdecomp ldiff_spe_dec = .false. ldiff_spe_dec_sp = .false. $end
```

应调用 `evaluatepara(init)`。

2. 两电极不同 scheme：

```fortran
$diffdecomp ldiff_spe_dec = .true. ldiff_spe_dec_sp = .false. $end
1   Pade   20
2   fano   3  3
```

应输出逐电极 scheme 信息，并调用 `evaluatepara_diff(init)`。

3. 大小写测试：

```text
Pade, pade, PSDJob, psdjob, Fano, fano, Prony, prony, MTSJOB, mtsjob, lprony, read
```

4. 非法输入：

- `ialf` 超出 `1:nalf`。
- 某个电极缺失配置。
- scheme 名称无法识别。
- pole 数小于等于 0。
- Fano/PSF 缺少正整数 type。

## 10. 当前限制

- diff 模式当前只支持 `ispectral=0` 的 Lorentzian/Drude 主路径。
- diff 模式暂不支持 `lband=.true.`。
- Fano 初始化仍复用原全局 `prelude_psdfff`，因此所有 Fano 电极必须使用同一个 `itype_fff`；若不同会报错停止。
- padding 维持旧数组结构，计算规模按最大 `ncor` 增大。
- `psdlor`、`pfdjob`、`mfdjob/nfreq` 未纳入本次逐电极混合支持。
