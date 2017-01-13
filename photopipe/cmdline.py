import json, sys, argparse, pkg_resources
from jsonschema import validate
from jsonschema.exceptions import ValidationError
from photopipe.reduction import preproc

def execute():
	# http://stackoverflow.com/a/23763283
	parser = argparse.ArgumentParser(prog='photopipe')
	parser.add_argument('-p', '--params', required=True, nargs='?', metavar='PARAMS.json', help='path of json file with reduction parameters')
	args = parser.parse_args()
#    collected_inputs = {'r': args.r}
	parameters_filename = args.params

	if parameters_filename:
		_parse_parameters_file(parameters_filename)

def _parse_parameters_file(filename):
	with open(filename, 'r') as params_file:		
		with open(pkg_resources.resource_filename(__name__, "schemas/parameters.json"), 'r') as schema_file:
			schema = json.load(schema_file)
			params = json.load(params_file)
			# print params
			# print schema
		try:
			validate(params, schema)
		except ValidationError as e:
			print e
			sys.exit(1)

		_execute_pipeline(params)


def _execute_pipeline(params):
	if 'preproc' in params:
		_execute_preprocessing(params['preproc'])

def _execute_preprocessing(params):
	# clean the params to sort out optionals
	if 'choose_calib' in params:
		choose_calib_calls = params['choose_calib']

		for call_params in choose_calib_calls:
			_execute_choose_calib(call_params)

def _execute_choose_calib(params):
	preproc.choose_calib( params['instrument'], 
						  params['ftype'],
						  params['workdir'],
						  params['cams'],
						  params['auto'],
						  params['reject_sat'],
						  params['amin'],
						  params['amax'],
						  params['save_select'],
						  params['figsize'],
						  params['noplot'] )



if __name__ == '__main__':
    execute()