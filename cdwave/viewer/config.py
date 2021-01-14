import json
import os
import logging
logger = logging.getLogger(__name__)
root = os.path.abspath(os.path.join(os.path.dirname(__file__)))

class Configure:
    

    def __init__(self):
        self.crt_config = {}
        self.all_config = self.load_config()
    
    def load_config(self):
        all_config = {'general': {}}
        if os.path.exists(os.path.join(root, '.config')):
            with open(os.path.join(root, '.config')) as fp:
                try:
                    all_config = json.load(fp)
                except ValueError:
                    logger.warning("Cannot parse the config file, need to be removed.")
                except ImportError:
                    pass
        return all_config

    def save_config(self):
        with open(os.path.join(root, '.config'), 'w') as fp:
            json.dump(self.all_config, fp, indent=2)

    def load(self, file_path, gui):
        if file_path in self.all_config:
            self.crt_config = self.all_config[file_path] 
        else:
            self.crt_config = {}
            self.all_config[file_path] = self.crt_config
        if 'filters' in self.crt_config:
            gui.load_filters(self.crt_config['filters'])

    def save(self, gui):
        self.all_config['general']['path'] = gui.lineEdit.text()
        if self.crt_config is not None:
            self.crt_config['filters'] = gui.filters
        self.save_config()
