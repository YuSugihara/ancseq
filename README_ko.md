# ancseq
#### 버전 1.2.1

## 목차
- [ancseq란 무엇인가?](#ancseq란-무엇인가)
- [설치](#설치)
  + [의존성](#의존성)
  + [conda를 사용한 설치](#conda를-사용한-설치)
- [사용법](#사용법)
  + [예제 1 : 뉴클레오타이드 서열 정렬을 위한 ancseq 실행](#예제-1--뉴클레오타이드-서열-정렬을-위한-ancseq-실행)
  + [예제 2 : 아미노산 서열 정렬을 위한 ancseq 실행](#예제-2--아미노산-서열-정렬을-위한-ancseq-실행)
  + [예제 3 : 코돈 서열 정렬을 위한 ancseq 실행](#예제-3--코돈-서열-정렬을-위한-ancseq-실행)
  + [예제 4 : 외부 그룹을 지정하여 ancseq 실행](#예제-4--외부-그룹을-지정하여-ancseq-실행)
  + [예제 5 : ```--fast``` 옵션으로 ancseq 실행](#예제-5---fast-옵션으로-ancseq-실행)
- [출력](#출력)
- [ancseq의 작업 흐름](#ancseq의-작업-흐름)
  + [IQ-TREE 명령어 1 : 계통수 구축](#iq-tree-명령어-1--계통수-구축)
  + [IQ-TREE 명령어 2 : 조상 서열 재구성](#iq-tree-명령어-2--조상-서열-재구성)
  + [IQ-TREE 명령어 3 : 조상 삽입 및 삭제 (INDEL) 재구성](#iq-tree-명령어-3--조상-삽입-및-삭제-indel-재구성)
- [anceseq는 DNA 모드에서 코돈 확률을 어떻게 계산하는가?](#anceseq는-dna-모드에서-코돈-확률을-어떻게-계산하는가)


## ancseq란 무엇인가?
<img src="https://github.com/YuSugihara/ancseq/blob/main/images/ancseq_workflow.png" width=400>

조상 서열 재구성은 다중 서열 정렬로부터 조상 상태를 재구성하는 기술입니다. ancseq는 [IQ-TREE](http://www.iqtree.org)를 사용하여 조상 서열을 재구성하는 래퍼 도구입니다. ancseq의 자세한 작업 흐름은 [여기](#ancseq의-작업-흐름)에서 확인하십시오.

#### 인용
Sugihara Y, Kourelis J, Contreras MP, Pai H, Harant A. Selvaraj M, Toghani A, Martinez-Anaya C*, Kamoun S* (2025) [Helper NLR immune protein NRC3 evolved to evade inhibition by a cyst nematode virulence effector](https://doi.org/10.1371/journal.pgen.1011653). _PLOS Genetics_, 21:e1011653 *Corresponding authors


## 설치
### 의존성
#### 소프트웨어
- [IQ-TREE](http://www.iqtree.org)

#### 파이썬 (>=3.5) 라이브러리
- [biopython](https://biopython.org)

### conda를 사용한 설치
[anaconda](https://www.anaconda.com)를 사용하여 의존성과 함께 ancseq를 설치할 수 있습니다.
```bash
git clone https://github.com/YuSugihara/ancseq.git
cd ancseq
conda env create -f ancseq.yml
conda activate ancseq
pip install .
```

## 사용법

```bash
$ ancseq -h
사용법: ancseq -s <ALIGNED_FASTA> -m <MODE> -o <OUT_DIR> [-t <INT>]

ancseq 버전 1.2.1

옵션:
  -h, --help         도움말 메시지를 표시하고 종료합니다.
  -s , --seq         FASTA 형식의 서열 정렬.
  -m , --mode        서열 유형. [DNA/AA/CODON]
  -o , --out         출력 디렉토리. 지정된 이름은 존재하지 않아야 합니다.
  -t , --threads     스레드 수. [4]
  -b , --bootstrap   부트스트랩 복제본 수. [1000]
  --max-report       동일한 위치에서 보고할 모호한 사이트의 최대 수. [5]
  --min-prob         모호한 사이트로 보고될 최소 확률. [0.05]
  --min-gap-prob     조상 상태를 갭으로 대체할 최소 확률. [0.5]
  --fast             IQ-TREE에서 -fast 옵션 사용 [FLASE]
  --model            IQ-TREE의 치환 모델 지정. IQ-TREE는 기본적으로 ModelFinder를 사용하여 최적의 치환 모델을 검색합니다 [MFP]
  --outgroup         IQ-TREE의 외부 그룹 지정. [None]
  --stop-codon-prob  DNA 모드에서 코돈 확률 계산 중지 [FLASE]
  --asr-only         트리 구축을 건너뛰고 조상 상태만 재구성 [FLASE]
  -v, --version      프로그램 버전 번호를 표시하고 종료합니다.
```

**노드의 조상 상태에 대한 오해를 피하기 위해 외부 그룹을 지정하는 것이 좋습니다. 자세한 내용은 [여기](#예제-4--외부-그룹을-지정하여-ancseq-실행)를 참조하십시오.**

+ [예제 1 : 뉴클레오타이드 서열 정렬을 위한 ancseq 실행](#예제-1--뉴클레오타이드-서열-정렬을-위한-ancseq-실행)
+ [예제 2 : 아미노산 서열 정렬을 위한 ancseq 실행](#예제-2--아미노산-서열-정렬을-위한-ancseq-실행)
+ [예제 3 : 코돈 서열 정렬을 위한 ancseq 실행](#예제-3--코돈-서열-정렬을-위한-ancseq-실행)
+ [예제 4 : 외부 그룹을 지정하여 ancseq 실행](#예제-4--외부-그룹을-지정하여-ancseq-실행)
+ [예제 5 : ```--fast``` 옵션으로 ancseq 실행](#예제-5---fast-옵션으로-ancseq-실행)


### 예제 1 : 뉴클레오타이드 서열 정렬을 위한 ancseq 실행
```
ancseq -s test_nuc.fasta \
       -m DNA \
       -o out_dir
```

`-s` : fasta 형식의 뉴클레오타이드 서열 정렬.

`-m` : 서열 유형.

`-o` : 출력 디렉토리의 이름. 지정된 이름은 존재하지 않아야 합니다.

### 예제 2 : 아미노산 서열 정렬을 위한 ancseq 실행

```bash
ancseq -s test_nuc.fasta \
       -m AA \
       -o out_dir
```

`-s` : fasta 형식의 아미노산 서열 정렬.

`-m` : 서열 유형.

`-o` : 출력 디렉토리의 이름. 지정된 이름은 존재하지 않아야 합니다.

### 예제 3 : 코돈 서열 정렬을 위한 ancseq 실행

**!!!경고!!!** IQ-TREE는 코돈 치환 모델을 구현합니다. 그러나 입력한 정렬에 따라 계통수를 구축하는 데 너무 오래 걸릴 수 있습니다. 이 경우 DNA 모드에서 ancseq를 실행하는 것이 좋습니다. anceseq는 DNA 모드에서 각 코돈의 확률을 계산할 수 있습니다.

```bash
ancseq -s test_codon.fasta \
       -m CODON \
       -o out_dir
```

`-s` : fasta 형식의 코돈 서열 정렬.

`-m` : 서열 유형.

`-o` : 출력 디렉토리의 이름. 지정된 이름은 존재하지 않아야 합니다.

### 예제 4 : 외부 그룹을 지정하여 ancseq 실행

외부 그룹을 지정하지 않고도 조상 상태를 재구성할 수 있습니다. 그러나 트리를 시각화할 때 노드의 조상 상태가 잘못 해석될 수 있습니다. 따라서 노드의 조상 상태에 대한 오해를 피하기 위해 외부 그룹을 지정하는 것이 좋습니다. IQ-TREE는 기본적으로 뿌리 있는 트리를 뿌리 없는 트리로 변환합니다.

```
ancseq -s test_nuc.fasta \
       -m DNA \
       --outgroup seq_id \
       -o out_dir
```

`-s` : fasta 형식의 뉴클레오타이드 서열 정렬.

`-m` : 서열 유형.

`--outgroup` : 외부 그룹의 서열 ID.

`-o` : 출력 디렉토리의 이름. 지정된 이름은 존재하지 않아야 합니다.

### 예제 5 : ```--fast``` 옵션으로 ancseq 실행

```bash
ancseq -s test_nuc.fasta \
       -m DNA \
       -o out_dir \
       --fast
```

`-s` : fasta 형식의 뉴클레오타이드 서열 정렬.

`-m` : 서열 유형.

`-o` : 출력 디렉토리의 이름. 지정된 이름은 존재하지 않아야 합니다.

`--fast` : IQ-TREE에서 ```-fast``` 옵션을 사용합니다.


## 출력
`OUT_DIR` 내부는 다음과 같습니다.
```
├── 00_tree
│  ├── 00_iqtree.err
│  ├── 00_iqtree.out
│  ├── test_nuc.fasta
│  ├── test_nuc.fasta.bionj
│  ├── test_nuc.fasta.ckp.gz
│  ├── test_nuc.fasta.contree
│  ├── test_nuc.fasta.iqtree
│  ├── test_nuc.fasta.log
│  ├── test_nuc.fasta.mldist
│  ├── test_nuc.fasta.model.gz
│  ├── test_nuc.fasta.splits.nex
│  └── test_nuc.fasta.treefile
├── 10_asr
│  ├── 10_iqtree.err
│  ├── 10_iqtree.out
│  ├── test_nuc.fasta
│  ├── test_nuc.fasta.ckp.gz
│  ├── test_nuc.fasta.iqtree
│  ├── test_nuc.fasta.log
│  ├── test_nuc.fasta.state.gz
│  └── test_nuc.fasta.treefile
├── 20_indels
│  ├── 20_iqtree.err
│  ├── 20_iqtree.out
│  ├── test_nuc.fasta.binary
│  ├── test_nuc.fasta.binary.ckp.gz
│  ├── test_nuc.fasta.binary.iqtree
│  ├── test_nuc.fasta.binary.log
│  ├── test_nuc.fasta.binary.state.gz
│  └── test_nuc.fasta.binary.treefile
└── 30_result
   ├── ancestral_state_result.treefile
   ├── ancestral_state_result.fasta
   ├── ancestral_state_result_with_gap.fasta
   ├── ancestral_state_result.sort.tsv
   ├── ancestral_state_result.tsv.gz
   └── ancestral_state_result.codon_prob.tsv.gz
```
- IQ-TREE로 재구성된 계통수는 `00_tree`에서 찾을 수 있습니다.
- 조상 서열 재구성 결과는 `30_result`에서 찾을 수 있습니다.
  + `ancestral_state_result.treefile`: 노드 레이블이 있는 계통수.
  + `ancestral_state_result.fasta`: 갭이 없는 조상 서열의 FASTA 파일.
  + `ancestral_state_result_with_gap.fasta`: 갭이 있는 조상 서열의 FASTA 파일.
  + `ancestral_state_result.sort.tsv` : 조상 상태의 확률.

## ancseq의 작업 흐름

<img src="https://github.com/YuSugihara/ancseq/blob/main/images/ancseq_workflow.png" width=400>

- [IQ-TREE 명령어 1](#iq-tree-명령어-1) : 계통수 구축.
- [IQ-TREE 명령어 2](#iq-tree-명령어-2) : 조상 서열 재구성.
- [IQ-TREE 명령어 3](#iq-tree-명령어-3) : 삽입 및 삭제 (INDEL) 재구성.

### IQ-TREE 명령어 1

```bash
iqtree -s ${INPUT_FASTA} \
       -st ${SEQ_TYPE} \
       -T ${NUM_THREADS} \
       -B ${NUM_BOOTSTRAP} \
       -m MFP \
       1> /OUT_DIR/00_tree/00_iqtree.out \
       2> /OUT_DIR/00_tree/00_iqtree.err
```

`-s` : fasta 형식의 서열 정렬.

`-st` : 서열 유형.

`-T` : 스레드 수.

`-B` : 초고속 부트스트랩 복제본 수.

`-m MFP` : 확장된 모델 선택 후 트리 추론.

ancseq에서 ```--fast``` 옵션을 지정하면 IQ-TREE 명령어 1이 다음과 같이 변경됩니다.

```bash
iqtree -s ${INPUT_FASTA} \
       -st ${SEQ_TYPE} \
       -T ${NUM_THREADS} \
       --alrt ${NUM_BOOTSTRAP} \
       -m MFP \
       --fast \
       1> /OUT_DIR/00_tree/00_iqtree.out \
       2> /OUT_DIR/00_tree/00_iqtree.err
```

`-s` : fasta 형식의 서열 정렬.

`-st` : 서열 유형.

`-T` : 스레드 수.

`--alrt` : SH 근사 우도비 검정 복제본 수.

`-m MFP` : 확장된 모델 선택 후 트리 추론.

`--fast` : FastTree와 유사한 빠른 검색.

### IQ-TREE 명령어 2

```bash
iqtree -asr \
       -s ${INPUT_FASTA} \
       -te /OUT_DIR/00_tree/${INPUT_FASTA}.treefile \
       -st ${SEQ_TYPE} \
       -T ${NUM_THREADS} \
       -m ${MODEL} \
       -o ${OUTGROUP} \
       -keep_empty_seq \
       1> /OUT_DIR/10_asr/10_iqtree.out \
       2> /OUT_DIR/10_asr/10_iqtree.err
```

`-asr` : 경험적 베이즈에 의한 조상 상태 재구성.

`-s` : fasta 형식의 서열 정렬.

`-te` : 트리 파일.

`-st` : 서열 유형.

`-T` : 스레드 수.

`-m` : 모델 이름.

`-o` : `--outgroup`으로 외부 그룹을 지정한 경우 외부 그룹의 서열 ID.

`-keep_empty_seq` : 정렬에서 빈 서열 유지.


### IQ-TREE 명령어 3

```bash
iqtree -asr \
       -s ${INPUT_FASTA}.binary \
       -te /OUT_DIR/00_tree/${INPUT_FASTA}.treefile \
       -st BIN \
       -T ${NUM_THREADS} \
       -blfix \
       -m JC2 \
       -o ${OUTGROUP} \
       -keep_empty_seq \
       1> /OUT_DIR/20_indels/20_iqtree.out \
       2> /OUT_DIR/20_indels/20_iqtree.err
```

`-asr` : 경험적 베이즈에 의한 조상 상태 재구성.

`-s` : fasta 형식의 서열 정렬.

`-te` : 트리 파일.

`-st BIN` : 이진 서열 유형.

`-T` : 스레드 수.

`-blfix` : `-t` 또는 `-te`를 통해 전달된 트리의 분기 길이 고정.

`-m JC2` : Jukes-Cantor 유형 이진 모델.

`-o` : `--outgroup`으로 외부 그룹을 지정한 경우 외부 그룹의 서열 ID.

`-keep_empty_seq` : 정렬에서 빈 서열 유지.

#### 참고문헌

1. Aadland K, Pugh C, Kolaczkowski B. 2019. High-Throughput Reconstruction of Ancestral Protein Sequence, Structure, and Molecular Function. In: Sikosek T ed. Computational Methods in Protein Evolution. Methods in Molecular Biology. New York, NY: Springer, 135–170. DOI: [10.1007/978-1-4939-8736-8_8](https://doi.org/10.1007/978-1-4939-8736-8_8).

2. VanAntwerp J, Finneran P, Dolgikh B, Woldring D. 2022. Ancestral Sequence Reconstruction and Alternate Amino Acid States Guide Protein Library Design for Directed Evolution. In: Traxlmayr MW ed. Yeast Surface Display. Methods in Molecular Biology. New York, NY: Springer US, 75–86. DOI: [10.1007/978-1-0716-2285-8_4](https://doi.org/10.1007/978-1-0716-2285-8_4).

## anceseq는 DNA 모드에서 코돈 확률을 어떻게 계산하는가?

IQ-TREE는 코돈 치환 모델을 구현합니다. 그러나 입력한 정렬에 따라 계통수를 구축하는 데 너무 오래 걸릴 수 있습니다. 따라서 각 뉴클레오타이드의 확률을 단순히 곱하여 각 코돈의 확률을 계산하는 함수를 구현했습니다. 예를 들어, $j$번째 위치에서 메티오닌의 확률 $P_{j, M}$은 $j$번째 코돈의 첫 번째 뉴클레오타이드에서 A의 확률 $p_{j_{1}, A}$, 두 번째 뉴클레오타이드에서 T의 확률 $p_{j_{2}, T}$, 세 번째 뉴클레오타이드에서 G의 확률 $p_{j_{3}, G}$를 곱하여 다음과 같이 계산할 수 있습니다.

$P_{j, M} = p_{j_{1}, A} * p_{j_{2}, T} * p_{j_{3}, G}$
