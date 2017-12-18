import numpy as np
import random
import cPickle
import load_data
import mscentipede
import argparse
import warnings
import time

# ignore warnings with these expressions
warnings.filterwarnings('ignore', '.*overflow encountered.*',)
warnings.filterwarnings('ignore', '.*divide by zero.*',)
warnings.filterwarnings('ignore', '.*invalid value.*',)

def learn_model(options):

    # write algorithm parameters
    log_handle = open(options.log_file,'w')
    log_handle.write('model = %s\n'%options.model)
    log_handle.write('window size = %d\n'%options.window)
    log_handle.write('motif file: %s\n'%options.motif_file)
    log_handle.write('bam files: %s\n'%(','.join([file for file in options.bam_files])))
    log_handle.close()

    # load motif sites
    log = "loading motifs ... "
    handle = open(options.log_file,'a')
    handle.write(log)
    handle.close()
    print log

    motif_handle = load_data.ZipFile(options.motif_file)
    locations = motif_handle.read()
    motif_handle.close()
    if np.any([len(loc)<5 for loc in locations]):
        print "error: ensure all rows in motif instance file contain same number of columns"
        sys.exit(1)

    if len(locations)>options.batch:
        order = np.argsort(np.array([loc[4] for loc in locations]).astype('float'))
        if order.size>3*options.batch:
            order = order[-3*options.batch:]
            random.shuffle(order)
            order = order[:options.batch]
        else:
            order = order[-options.batch:]
        locations = [locations[o] for o in order]

    try:
        scores = np.array([loc[4:] for loc in locations]).astype('float')
    except ValueError:
        print "error: column 5 and higher should all be numeric values."
        sys.exit(1)

    handle = open(options.log_file,'a')
    handle.write("done\n")
    handle.write("num of motif sites = %d\n"%(scores.shape[0]))
    handle.close()
    print "num of motif sites = %d"%(scores.shape[0])

    # load read data within specified window size
    log = "loading read counts ... "
    handle = open(options.log_file,'a')
    handle.write(log)
    handle.close()
    print log

    bam_handles = [load_data.BamFile(bam_file, options.protocol) for bam_file in options.bam_files]
    count_data = np.array([bam_handle.get_read_counts(locations, width=options.window) \
        for bam_handle in bam_handles])
    ig = [handle.close() for handle in bam_handles]   
    total_counts = np.sum(count_data, 2).T
    counts = np.array([count.T for count in count_data]).T

    handle = open(options.log_file,'a')
    handle.write("done\n")
    handle.close()

    # specify background
    if options.model=='msCentipede':
        if options.protocol=='DNase_seq':
            background_counts = np.ones((1,2*options.window,1), dtype=float)
        elif options.protocol=='ATAC_seq':
            background_counts = np.ones((1,options.window,1), dtype=float)
    elif options.model in ['msCentipede_flexbgmean','msCentipede_flexbg']:
        log = "Loading naked-DNA read counts ... "
        handle = open(options.log_file,'a')
        handle.write(log)
        handle.close()
        print log

        bam_handle = load_data.BamFile(options.bam_file_genomicdna, options.protocol)
        bg_count_data = np.array([bam_handle.get_read_counts(locations, width=options.window)])
        bam_handle.close()
        background_counts = np.array([count.T for count in bg_count_data]).T

        handle = open(options.log_file,'a')
        handle.write("done\n")
        handle.close()

    # estimate model parameters
    footprint_model, count_model, prior = mscentipede.estimate_optimal_model(counts, total_counts, scores, \
        background_counts, options.model, options.log_file, options.restarts, options.mintol)

    log = "writing model to file ... "
    handle = open(options.log_file,'a')
    handle.write(log)
    handle.close()
    print log

    # save model parameter estimates
    model_handle = open(options.model_file, 'w')
    cPickle.Pickler(model_handle,protocol=2).dump(footprint_model)
    cPickle.Pickler(model_handle,protocol=2).dump(count_model)
    cPickle.Pickler(model_handle,protocol=2).dump(prior)
    model_handle.close()

    handle = open(options.log_file,'a')
    handle.write("done\n")
    handle.close()

def infer_binding(options):

    # load the model
    handle = open(options.model_file, "r")
    footprint_model = cPickle.load(handle)
    count_model = cPickle.load(handle)
    prior = cPickle.load(handle)
    handle.close()

    # check if specified window size matches window size in model parameters
    if options.protocol=="DNase_seq":
        if 2**footprint_model[0].J!=options.window*2:
            print "Window size in model (%d bp) different from specified window size (%d bp). Using size in model ... \n"%(2**footprint_model[0].J/2, options.window)
            options.window = 2**footprint_model[0].J/2
    elif options.protocol=="ATAC_seq":
        if 2**footprint_model[0].J!=options.window:
            print "Window size in model (%d bp) different from specified window size (%d bp). Using size in model ... \n"%(2**footprint_model[0].J, options.window)
            options.window = 2**footprint_model[0].J

    # load motifs
    motif_handle = load_data.ZipFile(options.motif_file)
    
    # open read data handles
    bam_handles = [load_data.BamFile(bam_file, options.protocol) for bam_file in options.bam_files]

    # open background data handles
    if options.model=='msCentipede':
        if options.protocol=='DNase_seq':
            background_counts = np.ones((1,2*options.window,1), dtype=float)
        elif options.protocol=='ATAC_seq':
            background_counts = np.ones((1,options.window,1), dtype=float)
    elif options.model in ['msCentipede_flexbgmean','msCentipede_flexbg']:
        bg_handle = load_data.BamFile(options.bam_file_genomicdna, options.protocol)

    # check number of motif sites
    pipe = load_data.subprocess.Popen("zcat %s | wc -l"%options.motif_file, \
        stdout=load_data.subprocess.PIPE, shell=True)
    Ns = int(pipe.communicate()[0].strip())
    loops = Ns/options.batch+1

    # open gzip file to save inference
    handle = load_data.gzip.open(options.posterior_file, "wb")
    header = ['Chr','Start','Stop','Strand','LogPosOdds','LogPriorOdds','MultLikeRatio','NegBinLikeRatio']
    handle.write('\t'.join(header)+'\n')

    for n in xrange(loops):
        starttime = time.time()
        locations = motif_handle.read(batch=options.batch)

        count_data = np.array([bam_handle.get_read_counts(locations, width=options.window) \
            for bam_handle in bam_handles])
        total_counts = np.sum(count_data, 2).T
        counts = np.array([count.T for count in count_data]).T

        scores = np.array([loc[4:] for loc in locations]).astype('float')

        # specify background
        if options.model in ['msCentipede_flexbgmean','msCentipede_flexbg']:
            bg_count_data = np.array([bg_handle.get_read_counts(locations, width=options.window)])
            background_counts = np.array([count.T for count in bg_count_data]).T

        posterior_log_odds, prior_log_odds, footprint_log_likelihood_ratio, \
            total_log_likelihood_ratio = mscentipede.infer_binding_posterior(counts, \
            total_counts, scores, background_counts, \
            footprint_model, count_model, prior, options.model)
        output = np.hstack((posterior_log_odds, prior_log_odds, \
            footprint_log_likelihood_ratio, total_log_likelihood_ratio))

        towrite = [loc[:4] for loc in locations]

        ignore = [loc.extend(['%.4f'%p for p in pos])
            for loc,pos in zip(towrite,output)]
        ignore = [handle.write('\t'.join(map(str,elem))+'\n') for elem in towrite]
        print len(locations), time.time()-starttime

    handle.close()
    if options.model in ['msCentipede_flexbgmean','msCentipede_flexbg']:
        bg_handle.close()
    ig = [handle.close() for handle in bam_handles]


def parse_args():

    parser = argparse.ArgumentParser(description="runs msCentipede, to "
        "infer transcription factor binding, given a set of motif instances and "
        "chromatin accessibility data")

    parser.add_argument("--task",
                        choices=("learn","infer"),
                        default="learn",
                        help="specify whether to learn model parameters "
                        " or infer factor binding (default: learn)")

    parser.add_argument("--protocol",
                        choices=("ATAC_seq","DNase_seq"),
                        default="DNase_seq",
                        help="specifies the chromatin accessibility protocol (default:DNase_seq)")

    parser.add_argument("--model", 
                        choices=("msCentipede", "msCentipede_flexbg", "msCentipede_flexbgmean"),
                        default="msCentipede",
                        help="models differ in how they capture background rate of enzyme cleavage (default:msCentipede)")

    parser.add_argument("--restarts", 
                        type=int, 
                        default=1, 
                        help="number of re-runs of the algorithm (default: 1)")

    parser.add_argument("--mintol", 
                        type=float, 
                        default=1e-6,
                        help="convergence criterion for change in per-site marginal likelihood (default: 1e-6)")

    parser.add_argument("--model_file", 
                        type=str, 
                        default=None, 
                        help="file name to store the model parameters")

    parser.add_argument("--posterior_file", 
                        type=str, 
                        default=None, 
                        help="file name to store the posterior odds ratio, and "
                        "likelihood ratios for each model component, at each motif. ")

    parser.add_argument("--log_file",
                        type=str,
                        default=None,
                        help="file name to store some statistics of the EM algorithm ")
#                        "and a plot of the cleavage profile at bound sites")

    parser.add_argument("--window", 
                        type=int, 
                        default=128, 
                        help="size of window around the motif instance, where chromatin "
                        "accessibility profile is used for inferring transcription "
                        "factor binding. (default: 128)")

    parser.add_argument("--batch", 
                        type=int, 
                        default=10000, 
                        help="maximum number of motif instances used for learning model parameters. "
                        " this is also the number of motif instances on which inference is "
                        " performed at a time. (default: 10000)")

    parser.add_argument("motif_file",
                        action="store",
                        help="name of a gzipped text file containing "
                        " positional information and other attributes for motif instances "
                        " of a transcription factor. columns of the file should be as follows. "
                        " Chromosome Start End Strand PWM_Score [Attribute_1 Attribute_2 ...]. "
                        " additional attributes are optional.")

    parser.add_argument("bam_files",
                        action="store",
                        nargs="+",
                        help="whitespace separated list of bam files "
                        " from a chromatin accessibility assay ")

    parser.add_argument("--bam_file_genomicdna",
                        action="store",
                        default=None,
                        help="bam file from a chromatin accessibility assay on genomic DNA")

    parser.add_argument("--seed",
                        default=None,
                        help="set seed for random initialization of parameters")

    options = parser.parse_args()

    # if no motif file is provided, throw an error
    if options.motif_file is None:
        parser.error("Need to provide a file of motifs for a transcription factor")

    # if no model file is provided, create a `default` model file name
    if options.model_file is None:
        options.model_file = "%s_%s_model_parameters.pkl"%(options.motif_file.split('.')[0], '_'.join(options.model.split('-')))

    # if no posterior file is provided, create a `default` posterior file name
    if options.posterior_file is None:
        options.posterior_file = "%s_%s_binding_posterior.txt.gz"%(options.motif_file.split('.')[0], '_'.join(options.model.split('-')))

    # if no log file is provided, create a `default` log file name
    if options.log_file is None:
        options.log_file = "%s_%s_log.txt"%(options.motif_file.split('.')[0],'_'.join(options.model.split('-')))
    
    # make sure model file exists, before trying to run inference
    if options.task=='infer':
        try:
            handle = open(options.model_file, 'r')
            handle.close()
        except IOError:
            parser.error("Need to provide the file where model parameters are saved")

    if options.model in ['msCentipede_flexbgmean','msCentipede_flexbg'] and options.bam_file_genomicdna is None:
        parser.error("Need to provide a bam file containing chromatin accessibility "
            "data in genomic DNA, if the model is specified to be "
            "msCentipede-flexbgmean or msCentipede-flexbg")

    if options.seed is not None:
        np.random.seed(int(options.seed))

    if options.window<=0:
        options.window = 128

    if not ((options.window & (options.window - 1)) == 0):
        options.window = 2**(int(np.log2(options.window)))

    return options


def main():

    options = parse_args()

    if options.task=='learn':
        learn_model(options)

    elif options.task=='infer':
        infer_binding(options)

if __name__=="__main__":

    main()
