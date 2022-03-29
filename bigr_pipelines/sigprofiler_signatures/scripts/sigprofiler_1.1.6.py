
import sys
import os
import os.path


# /mnt/beegfs/software/conda/bin/pip install --user matplotlib==3.4.3
# /mnt/beegfs/software/conda/bin/pip install --user SigProfilerMatrixGenerator

# installation of references in /home/m_diop@intra.igr.fr/.local/lib/python3.7/site-packages/SigProfilerMatrixGenerator/references
# /mnt/beegfs/software/conda/bin/python -c \
if "--install" in sys.argv:
    from SigProfilerMatrixGenerator import install as genInstall
    genInstall.install(sys.argv[1], bash=True)

# /mnt/beegfs/software/conda/bin/pip install --user SigProfilerExtractor
if os.path.exists(sys.argv[2]):
    os.chdir(sys.argv[2])
    print("Working in {}".format(os.getcwd()))
    print(os.listdir("vcf"))

#if not os.path.exists("vcf/logs/"):
#    os.makedirs("vcf", exist_ok = True)
#    os.makedirs("vcf/logs", exist_ok = True)
#    print("vcf/logs/ had to be created")
    
#if not os.path.exists("sigprofiler/vcf/logs/"):
#    os.makedirs("sigprofiler", exist_ok = True)
#    os.makedirs("sigprofiler/vcf", exist_ok = True)
#    os.makedirs("sigprofiler/vcf/logs", exist_ok = True)
#    print("sigprofiler/vcf/logs had to be created")
    


from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor(
    'vcf',
    'results',
    'vcf',
    reference_genome=sys.argv[1],
    exome=True,
    minimum_signatures=1,
    maximum_signatures=10,
    nmf_replicates=100,
    cpu=10
)
