import json
import airr

def create_json_job_files(data, project_id, save=True):
    '''
    Parameters
    ----------
    data : dict
        - From `airr.read_airr('/path/to/repertoires.airr.json')`
    project_id : str
        - The project's UUID
    save : boolean
        -  If `True`, JSON will be dumped and make one job file per library.

    Returns
    -------
    json_contents : list
        - Returns a list where each element is JSON job file contents for a library from `data`
    '''
    libraries = []
    curr_l_num = -1
    for rep in data['Repertoire']:
        if curr_l_num != rep['sample'][0]['sequencing_run_id']:
            libraries.append([])
            curr_l_num = rep['sample'][0]['sequencing_run_id']
        libraries[-1].append((rep['sample'][0]['sequencing_run_id'], rep['sample'][0]['sequencing_files']))

    seq_for_files_source_urls = []
    seq_for_files = []
    seq_rev_files_source_urls = []
    seq_rev_files = []
    for lib in libraries:
        seq_for_files_source_urls.append([])
        seq_for_files.append([])
        seq_rev_files_source_urls.append([])
        seq_rev_files.append([])
        for i in range(len(lib)):
            seq_for_files_source_urls[-1].append('tapis://data-storage.vdjserver.org/projects/'+project_id+'/files/'+lib[i][1]['filename'])
            seq_for_files[-1].append(lib[i][1]['filename'])
            seq_rev_files_source_urls[-1].append('tapis://data-storage.vdjserver.org/projects/'+project_id+'/files/'+lib[i][1]['paired_filename'])
            seq_rev_files[-1].append(lib[i][1]['paired_filename'])
            
    seq_for_files = [' '.join(seq_for_files[i]) for i in range(len(seq_for_files))]
    seq_rev_files = [' '.join(seq_rev_files[i]) for i in range(len(seq_rev_files))]

    json_contents = []
    for i in range(len(libraries)):
        lib_num = libraries[i][0][0]
        content = {
            'name' : 'COVVAX - library ' + lib_num,
            'appId' : 'vdjpipe-ls6',
            'appVersion' : '0.2',
            'maxMinutes' : 8*60,
            'nodeCount' : 1,
            'archiveSystemId' : 'data-storage.vdjserver.org',
            'archiveSystemDir' : '/projects/'+project_id+'/analyses/${JobUUID}',
            'fileInputArrays' : [{'name' : 'SequenceForwardPairedFiles', 'sourceUrls' : seq_for_files_source_urls[i]},
                                {'name' : 'SequenceReversePairedFiles', 'sourceUrls' : seq_rev_files_source_urls[i]}],
            'fileInputs' : [{'name' : 'ForwardPrimerFile', 'sourceUrl' : 'tapis://data-storage.vdjserver.org/projects/'+project_id+'/files/primers.fasta', 'targetPath' : 'primers.fasta'}],
            'parameterSet' : {
                'schedulerOptions' : [
                    {'name' : 'allocation', 'arg' : '-A MCB23006'}
                ],
                'containerArgs' : [
                ],
                'appArgs' : [
                ],
                'envVariables' : [
                    {'key' : 'Workflow', 'value' : 'paired'},
                    {'key' : 'SequenceForwardPairedFiles', 'value' : seq_for_files[i]},
                    {'key' : 'SequenceForwardPairedFilesMetadata', 'value' : '12348'},
                    {'key' : 'SequenceReversePairedFiles', 'value' : seq_rev_files[i]},
                    {'key' : 'SequenceReversePairedFilesMetadata', 'value' : '4567'},
                    {'key' : 'MergeMinimumScore', 'value' : '10'},
                    {'key' : 'PreFilterStatisticsFlag', 'value' : '1'},
                    {'key' : 'FilterFlag', 'value' : '1'},
                    {'key' : 'PostFilterStatisticsFlag', 'value' : '1'},
                    {'key' : 'MinimumAverageQuality', 'value' : '35'},
                    {'key' : 'MinimumLength', 'value' : '200'},
                    {'key' : 'MaximumHomopolymer', 'value' : '20'},
                    {'key' : 'ForwardPrimer', 'value': '1'},
                    {'key' : 'ForwardPrimerMaximumMismatches', 'value': '1'},
                    {'key' : 'ForwardPrimerTrim', 'value': '1'},
                    {'key' : 'ForwardPrimerSearchWindow', 'value': '50'},
                    {'key' : 'FindUniqueFlag', 'value' : '1'}
                ]
            }
        }
        json_contents.append(content)

    file_names = ['job-vdjpipe-COVVAX-library-'+lib[0][0]+'.json' for lib in libraries]
    
    if save:
        for i in range(len(file_names)):
            with open(file_names[i], 'w') as json_file:
                json.dump(json_contents[i], json_file, indent=4)
    
    return json_contents

def main():
    project_id='c0aef9ef-362d-46fd-abba-39ac0945bd66'
    data = airr.read_airr('repertoires_COVVAX.airr.json')

    create_json_job_files(data, project_id)

if __name__ == '__main__':
    main()