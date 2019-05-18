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
from Bio import motifs,SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord


app = Flask(__name__)

# model = load_model('static/repeat_validator_model.h5')
# model._make_predict_function()
# graph = tf.get_default_graph()
model = load_model('static/repeat_validator_43_new_ds_16_e.h5')
input_len = 43
validation_threshold = 0.7

# loading
with open('static/dna_tokenizer.pickle', 'rb') as handle:
    tok = pickle.load(handle)

def randomDigits(digits):
    lower = 10**(digits-1)
    upper = 10**digits - 1
    return random.randint(lower, upper)

def get_validation_score(repeat):
   tsp = tok.texts_to_sequences([repeat])
   tsp_matrix = sequence.pad_sequences(tsp,maxlen=input_len)
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
      if request.form.get('fasta-sequence')!='':
         filename = str(randomDigits(9))+".fasta"
         with open('uploaded_sequences/'+filename,"w+") as f:
            f.write(request.form.get('fasta-sequence'))
      else:
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
   try:
      with open('output'+file_name+'.json') as f:
         data = json.load(f)
   except:
      return jsonify({'error':'No CRISPR arrays found in the input sequence!'})
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
def create_webLogo(file_name):
    repeats = request.get_json()
    m = motifs.create(list(repeats))
    m.weblogo('static/logos/'+str(file_name)+'.png')
    return jsonify('{"success":1}')

@app.route('/generate_dna_struct/<file_name>', methods = ['POST'])
def create_dna_structure(file_name):
   results = request.get_json()
   features = []
   for i, spacerRepeat in enumerate(results['spacerRepeats']):
      features.append(GraphicFeature(start=spacerRepeat['position'], end=spacerRepeat['position']+len(spacerRepeat['repeat']), 
                  strand=+1, color="#cffccc", label="Repeat_"+str(i+1)))
      if 'spacer' in spacerRepeat:
         features.append(GraphicFeature(start=spacerRepeat['position']+len(spacerRepeat['repeat'])+1, 
         end=spacerRepeat['position']+len(spacerRepeat['repeat'])+spacerRepeat['lengths'][1], strand=+1, color="#ccccff",
                   label="Spacer_"+str(i+1)))
   record = GraphicRecord(sequence_length=results['length'], features=features)
   record = record.crop((results['spacerRepeats'][0]['position']-50, 
      results['spacerRepeats'][len(results['spacerRepeats'])-1]['position']+
      len(results['spacerRepeats'][len(results['spacerRepeats'])-1]['repeat'])+50))
   ax, _ = record.plot(figure_width=10)
   ax.figure.savefig('static/logos/'+str(file_name)+'.png', bbox_inches='tight')
   return jsonify('{"success":1}')
		
@app.route('/generate_structure/<file_name>', methods = ['POST'])
def create_2dStructure(file_name):
   sequence = request.get_json()
   with open("uploaded_sequences/"+file_name+".fasta", "w") as output_handle:
      output_handle.write(">repeat\n"+sequence['seq'])

   if os.path.isfile('uploaded_sequences/'+file_name+'.fasta'):
      print('file found', 'uploaded_sequences/'+file_name+'.fasta')
      subprocess.call(['C:\Program Files\RNAstructure6.1\exe\Fold', 'uploaded_sequences/'+file_name+'.fasta', 'uploaded_sequences/'+file_name+'.ct',
       '--DNA', '--loop' , '30', '--maximum', '20', '--percent', '10', '--temperature', '310.15', '--window', '3'])
   else:
      raise Exception()

   if os.path.isfile('uploaded_sequences/'+file_name+'.ct'):
      subprocess.call(['C:\Program Files\RNAstructure6.1\exe\draw', 'uploaded_sequences/'+file_name+'.ct', 'static/logos/'+file_name+'.svg',
       '--svg', '-n' , '1'])
   else:
      raise Exception()

   return jsonify('{"success":1}')

@app.route('/get_content/<file_name>', methods = ['GET'])
def get_content(file_name):
   content = ''
   try:
      with open("static/examples/"+file_name+".fasta", 'r') as f:
         content = f.read()
   except Exception as e:
      return 'Getting content error: '+str(e),500
   return content

if __name__ == '__main__':
   app.run(host='0.0.0.0', debug=False, port=80)
  