% Imputation script
% Ref

%% Requrement file
    % Reference data
    mkdir 01_rawdata
    copyfile('../00_rawdata/v3.*', '01_rawdata');
    cd('01_rawdata');
    system('tar -xzvf v3.20101123.ENIGMA2.EUR.20120719.vcf.tgz');
    system('tar -xzvf v3.20101123.ENIGMA2.EUR.20120719.extras.tgz');

    % Input data
    input_file_name = 'ppmi_genotypes_filt_CEU';
    copyfile(['../../00_rawdata/' input_file_name '*'],'../01_rawdata');

%% Preprocessing
    cd ..
    mkdir 02_preprocess
    cd 02_preprocess
    
    copyfile(['../01_rawdata/' input_file_name '*'],'../02_preprocess');

    % re-quailtity control
    re_QC_command_1 = ['awk ''{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}'' ' input_file_name '.bim | grep ambig | awk ''{print $1}'' > ambig.list'];
    system(re_QC_command_1);
    
    re_QC_command_2 = ['plink --bfile ' input_file_name ' --exclude ambig.list --make-founders --out lastQC --maf 0.01 --hwe 0.000001 --make-bed --noweb'];
    system(re_QC_command_2);
    
    % flipping
    reformat_command_1 = 'awk ''{print $2,$1,$3,$4,$5,$6}'' lastQC.bim > tempQC.bim';
    reformat_command_2 = 'awk ''NR==FNR{s=$1;a[s]=$0;next} a[$1]{print $0 " "a[$1]}'' tempQC.bim ../01_rawdata/1kgp.alleles > merged.alleles';    
    reformat_command_3 = 'awk ''{ if ($2!=$8 && $2!=$9) print $1}'' merged.alleles > flip.list';       
    reformat_command_4 = 'plink --bfile lastQC --extract ../01_rawdata/1kgp.snps --update-map ../01_rawdata/1kgp.chr --update-chr --flip flip.list --make-bed --out temp --noweb';        
    reformat_command_5 = 'plink --bfile temp --update-map ../01_rawdata/1kgp.bp --geno 0.05 --mind 0.05 --make-bed --out lastQCb37 --noweb';
    reformat_command_6 = 'wc -l lastQCb37.bim';
    system(reformat_command_1); system(reformat_command_2); system(reformat_command_3); system(reformat_command_4); system(reformat_command_5); system(reformat_command_6);
    
    % X-chromosome
    reformat_command_7 = 'awk ''{ if($5==1) print $1, $2}'' lastQCb37.fam > male.list';
    reformat_command_8 = 'awk ''{ if($5==2) print $1, $2}'' lastQCb37.fam > female.list';
    system(reformat_command_7); system(reformat_command_8);
    
%% MACH
    % split chrmosomes
    disp('Split chromosom');
    mkdir ../03_Mach
    cd ../03_Mach
    mkdir 00_split_chr
    copyfile('../00_data/00_split_chr.sh','../03_Mach');
    system('chmod +x 00_split_chr.sh');
    system('./00_split_chr.sh');
    movefile('00_split_chr/ready4mach.23.female.map', '00_split_chr/ready4mach.23.map');
    delete('00_split_chr/ready4mach.23.male.map');
    
    % reformat the map file to make a merlin dat file in which the SNP are named as chromosome:position format
    disp('Reformatting data');
    mkdir 01_mach_reformat
    copyfile('../00_data/01_mach_reformat.sh','../03_Mach');
    system('chmod +x 01_mach_reformat.sh');
    system('./01_mach_reformat.sh'); 
    copyfile('00_split_chr/ready4mach.*.ped','01_mach_reformat');
    cd 01_mach_reformat
    system('gzip ready4mach.*.ped');
    
    % MACH process
    mkdir ../02_Mach_process
    copyfile('../../00_data/02_MaCH_phasing_1.sh','./');
    system('chmod +x 02_MaCH_phasing_1.sh');
    system('./02_MaCH_phasing_1.sh');
    movefile('autoChunk-*','../02_Mach_process');
    movefile('chunk*','../02_Mach_process');
    
    cd ..
    mkdir 02_Mach_process/tmp
    copyfile('../00_data/02_MaCH_phasing_2.sh','./');
    system('chmod +x 02_MaCH_phasing_2.sh');
    system('./02_MaCH_phasing_2.sh');
    
    % MACH process - execute
    copyfile('../00_data/03_MaCH_job.txt','./');
    system('chmod +x 03_MaCH_job.txt');
    system('./03_MaCH_job.txt');
    
    % MACH for x chromosome
    x_mach = 'zless 01_mach_reformat/ready4mach.23.male.ped.gz | awk ''{ printf "%s", $1 "->"$2 " HAPLO1 "; for(N=7; N<=NF; N+=2) printf "%s", $N; printf "\n"; printf "%s", $1 "->"$2 " HAPLO2 "; for(N=7; N<=NF; N+=2) printf "%s", $N; printf "\n";}'' > 02_Mach_process/ready4mach.23.male';
    system(x_mach);
    system('gzip -f 02_Mach_process/ready4mach.23.male > 02_Mach_process/ready4mach.23.male.gz');
    
 %% Minimac
    % MiniMac process
    mkdir 03_Minimac_process
    mkdir 03_Minimac_process/tmp
    copyfile('../00_data/04_MiniMac-impute.sh','./');
    system('chmod +x 04_MiniMac-impute.sh');
    system('./04_MiniMac-impute.sh');
    
    % MiniMac process - execute
    copyfile('../00_data/05_Minimac_job.txt','./');
    system('chmod +x 05_Minimac_job.txt');
    system('./05_Minimac_job.txt');
    
 %% Merge and re-format
    % Convert format
    mkdir 04_merged_imputation
    mkdir 04_merged_imputation/tmp
    copyfile('../00_data/06_dos2plink.sh','./');
    system('chmod +x 06_dos2plink.sh');
    system('./06_dos2plink.sh');
    
    copyfile('../00_data/07_job_dose2plink.txt','./');
    system('chmod +x 07_job_dose2plink.txt');
    system('./07_job_dose2plink.txt')
    
    % Merge file
    merged_file_name = 'adni.imputed'; %% Change if necesseary
    cd 04_merged_imputation
    copyfile('../../00_data/list.txt','./');
    merge_command = ['plink --bfile chunk1.chr1.imputed --make-bed --merge-list list.txt --out ' merged_file_name ' --noweb'];
    system(merge_command);
    
    % Reformate data ->> user parameter 
    % 1) delete ambig snp
    % 2) reformat chr:bp > rs id
    reformat_output = 'reformat_output';
    re_qc = ['awk ''{ if ($5=="T"||$5=="A"||$5=="C"||$5=="G") print $2 ; else print $2, "ambig" ;}'' ' merged_file_name '.bim | grep ambig | awk ''{print $1}'' > ambig.list'];
    reformat = ['plink --bfile ' merged_file_name ' --exclude ambig.list --make-founders --out reformat_1 --maf 0.01 --hwe 0.000001 --make-bed --noweb'];
    system(re_qc);system(reformat);    
    
    copyfile('../../00_data/reform.py','./');
    system('chmod +x reform.py');
    system('python reform.py')
  
    copyfile('reformat_1.bed','reformat.bed');
    copyfile('reformat_1.fam','reformat.fam');
    
    mkdir ../final_imputation
    copyfile('../../00_data/rsID.txt','./');
    system('plink --bfile reformat --make-bed --out ../final_imputation/adni.imputed --update-map rsID.txt --update-name --noweb')
    
    
    
    
    
    
    
    
    
    
    
    
    