import os
os.environ['KERAS_BACKEND'] = 'theano'
from flask import Flask, render_template, request
from werkzeug import secure_filename
import time
import random
import subprocess
import os.path
from keras.models import load_model
import pickle
from keras.preprocessing import sequence
import json
from flask import jsonify
from Bio.Seq import Seq
from Bio import motifs

# import tensorflow as tf

app = Flask(__name__)

# model = load_model('static/repeat_validator_model.h5')
# model._make_predict_function()
# graph = tf.get_default_graph()
model = load_model('static/repeat_validator_model.h5')

# loading
with open('static/dna_tokenizer.pickle', 'rb') as handle:
    tok = pickle.load(handle)
validation_threshold = 0.7

def randomDigits(digits):
    lower = 10**(digits-1)
    upper = 10**digits - 1
    return random.randint(lower, upper)

def get_validation_score(repeat):
   tsp = tok.texts_to_sequences([repeat])
   tsp_matrix = sequence.pad_sequences(tsp,maxlen=45)
   # with graph.as_default():
      # score = model.predict(tsp_matrix)
   score = model.predict(tsp_matrix)
   return score[0][0]

'''
Referred from https://www.tutorialspoint.com/flask/flask_file_uploading.htm
'''
@app.route('/')
def render_index():
   return render_template('index.html')
	
@app.route('/uploader', methods = ['POST'])
def upload_file():
   if request.method == 'POST':
      try:
         f = request.files['seqfile']
         filename = str(randomDigits(9))+"."+f.filename.split('.')[-1]
         f.save('uploaded_sequences/'+secure_filename(filename))
      except:
         return 'File Upload Error', 400
      return filename.split('.')[0]

@app.route('/get_candidates/<file_name>', methods = ['GET'])
def get_candidates(file_name):
   try:
      if os.path.isfile('uploaded_sequences/'+file_name+'.fasta'):
         subprocess.call(['java', '-jar', 'modifiedCRT.jar', '-requestID', str(int(file_name)) ,'uploaded_sequences/'+file_name+'.fasta'])
      else:
         raise Exception()
   except Exception as e:
      return 'CRT Candidate Extraction error'+str(e),500
   return file_name

@app.route('/filter_candidates/<file_name>', methods = ['GET'])
def filter_candidates(file_name):
   with open('output'+file_name+'.json') as f:
      data = json.load(f)
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
   return jsonify(data)

@app.route('/generate_image/<file_name>', methods = ['POST'])
def your_route(file_name):
    repeats = request.get_json()
    m = motifs.create(list(repeats))
    m.weblogo('static/logos/'+str(file_name)+'.png')
    return jsonify('{"success":1}')
		
if __name__ == '__main__':
   app.run(debug = True)