import os
os.environ['KERAS_BACKEND'] = 'theano'
import time
import subprocess
import os.path
from keras.models import load_model
import pickle
from keras.preprocessing import sequence
import json
from Bio.Seq import Seq
from Bio import motifs,SeqIO
import time
import sys
from tqdm import tqdm

#load model
model = load_model('static/repeat_validator_model.h5')

# load tokenizer
with open('static/dna_tokenizer.pickle', 'rb') as handle:
    tok = pickle.load(handle)

# utility function to create output file
def create_output(data, fileName):
    with open(fileName, 'w+') as f:
        f.write(data)

# utility function to score repeat
def get_validation_score(repeat):
    tsp = tok.texts_to_sequences([repeat])
    tsp_matrix = sequence.pad_sequences(tsp, maxlen=45)
    score = model.predict(tsp_matrix)
    return score[0][0]


def runCRISPRLstm(input_fasta, output_txt, validation_threshold):
    
    if not os.path.isfile(input_fasta):
        print("Error, input fasta file not found!!")
        exit()

    # Create CRISPR candidates using linent CRT
    subprocess.call(['java', '-jar', 'modifiedCRT.jar', '-requestID', "100" , input_fasta])

    #load the crt candidate output json
    try:
        with open('output100.json') as f:
            data = json.load(f)
    except:
        create_output({'error':'No CRISPR arrays found in the input sequence!'}, output_txt)
        exit()

    validCrisprs=0
    start_time = time.time()
    for crispr in data['CRISPRs']:
        total_score = 0
        for spacer_repeat in crispr['spacerRepeat']:
            total_score += get_validation_score(spacer_repeat['repeat'])
        total_score/=crispr['noOfRepeats']
        if total_score >= validation_threshold:
            instances = []
            crispr['isValid'] = True
            validCrisprs+=1
            for spacer_repeat in crispr['spacerRepeat']:
                instances.append(spacer_repeat['repeat'])
            m = motifs.create(instances)
            pwm = m.counts.normalize(pseudocounts=0.5)
            crispr['consensusRepeat'] = str(pwm.consensus)
        else:
            crispr['isValid'] = False
    end_time = time.time()
    data['validCrisprs'] = validCrisprs
    data['timeTaken']+=end_time-start_time

    create_output(json.dumps(data), output_txt)
    os.remove('output100.json')



if __name__ == '__main__':  
    if len(sys.argv)<2:
        print("Error, please specificy input and output file arguments!")
        exit()

    input_fasta = sys.argv[1]
    output_txt = sys.argv[2]
    parameters = {
        "thresh":0.5,
        "dir":False
    }
    
    if len(sys.argv) > 3:
        i=3
        while i<len(sys.argv):
            parameters[sys.argv[i]] = sys.argv[i+1]
            i+=2
    
    print(parameters)
    if parameters["dir"]:
        for fasta_file in tqdm(os.listdir(input_fasta)):
            if '.fasta' not in fasta_file:
                continue
            runCRISPRLstm(input_fasta+"/"+fasta_file, output_txt+"/"+fasta_file.replace(".fasta",".txt"), parameters["thresh"])
    else:
        runCRISPRLstm(input_fasta, output_txt, parameters["thresh"])

"""
to run use command:
Bare minimum: python cli.py input_fasta_file output_txt_file
Specify threshold: python cli.py input_fasta_file output_txt_file thresh 0.5 
Specify threshold & directory mode: python cli.py dir_with_fastas output_dir thresh 0.5 dir True
"""
        
