"""Reusable JSON configuration storage and setup dialogs for iolite exporters."""

from copy import deepcopy
from dataclasses import dataclass
import json
import math
import os
import stat
import tempfile
from typing import Optional

from iolite.QtCore import QSettings
from iolite.QtGui import QFileDialog, QMessageBox


@dataclass(frozen=True)
class ConfigField:
    value_type: type
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    decimals: Optional[int] = None

    def validate(self, key, value):
        if self.value_type is list:
            if not isinstance(value, list) or not value or not all(
                isinstance(item, str) and item.strip() for item in value
            ):
                raise ValueError(f'"{key}" must be a non-empty list of non-empty strings.')
        elif self.value_type is str:
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f'"{key}" must be a non-empty string.')
        elif self.value_type is bool:
            if not isinstance(value, bool):
                raise ValueError(f'"{key}" must be a boolean.')
        elif self.value_type in (int, float):
            valid_type = (
                isinstance(value, int) if self.value_type is int
                else isinstance(value, (int, float))
            )
            if isinstance(value, bool) or not valid_type:
                raise ValueError(f'"{key}" must be a number of type {self.value_type.__name__}.')
            if self.minimum is not None and value < self.minimum:
                raise ValueError(f'"{key}" must be at least {self.minimum}.')
            if self.maximum is not None and value > self.maximum:
                raise ValueError(f'"{key}" must be at most {self.maximum}.')
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f'"{key}" must be finite.')
            if self.decimals is not None and round(value, self.decimals) != value:
                raise ValueError(f'"{key}" must have at most {self.decimals} decimal places.')
        else:
            raise ValueError(f'Unsupported field type for "{key}": {self.value_type}.')


class ExporterConfig:
    """Manage validated JSON configuration and setup dialogs for an exporter."""

    def __init__(self, defaults, fields, settings_key, title, choices=None):
        self.defaults = deepcopy(defaults)
        self.fields = dict(fields)
        self.settings_key = settings_key
        self.title = title
        self.choices = dict(choices or {})
        if self.defaults.keys() != self.fields.keys():
            raise ValueError('Defaults and fields must define the same configuration keys.')
        for selected_key, options_key in self.choices.items():
            if (
                selected_key not in self.fields or options_key not in self.fields
                or self.fields[selected_key].value_type is not str
                or self.fields[options_key].value_type is not list
            ):
                raise ValueError('Each choice must link a string field to a string-list field.')
        self.validate(self.defaults)

    def validate(self, overrides):
        if not isinstance(overrides, dict):
            raise ValueError('The configuration must be a JSON object.')
        unknown = overrides.keys() - self.fields.keys()
        if unknown:
            raise ValueError(f'Unknown configuration keys: {", ".join(sorted(unknown))}')
        config = deepcopy({**self.defaults, **overrides})
        for key, field in self.fields.items():
            field.validate(key, config[key])
        for selected_key, options_key in self.choices.items():
            if config[selected_key] not in config[options_key]:
                if selected_key in overrides:
                    raise ValueError(f'"{selected_key}" must be one of the values in "{options_key}".')
                config[selected_key] = config[options_key][0]
        return config

    def read_overrides(self, path):
        with open(path, encoding='utf-8') as config_file:
            overrides = json.load(config_file)
        self.validate(overrides)
        return overrides

    def load(self, path):
        return self.validate(self.read_overrides(path))

    def write(self, path, overrides):
        self.validate(overrides)
        path = os.path.realpath(path)
        contents = json.dumps(overrides, indent=4, ensure_ascii=False, allow_nan=False) + '\n'
        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode='w', encoding='utf-8', dir=os.path.dirname(path),
                prefix='.exporter-config-', suffix='.tmp', delete=False
            ) as config_file:
                temporary_path = config_file.name
                config_file.write(contents)
            if os.path.exists(path):
                os.chmod(temporary_path, stat.S_IMODE(os.stat(path).st_mode))
            os.replace(temporary_path, path)
        finally:
            if temporary_path is not None and os.path.exists(temporary_path):
                os.remove(temporary_path)

    def remember(self, path, parent=None):
        settings = QSettings()
        settings.setValue(self.settings_key, path)
        settings.sync()
        if settings.status() != QSettings.NoError:
            QMessageBox.warning(
                parent, 'Could not remember configuration',
                'The configuration is loaded, but its path could not be saved in '
                'iolite preferences. You may be asked to select it next time.'
            )

    def choose(self, parent=None, use_remembered=True):
        if use_remembered:
            path = QSettings().value(self.settings_key, '')
            if path:
                try:
                    return path, self.load(path)
                except (OSError, ValueError) as error:
                    QMessageBox.warning(
                        parent, 'Could not load configuration',
                        f'Could not load {path}:\n{error}\n\nSelect another file, create one, or cancel.'
                    )

        while True:
            message = QMessageBox(parent)
            message.setWindowTitle(f'{self.title} Configuration')
            message.setText('Select an existing configuration or create a new one using the shipped defaults.')
            select_button = message.addButton('Select config file', QMessageBox.ActionRole)
            create_button = message.addButton('Create new config file', QMessageBox.ActionRole)
            message.addButton(QMessageBox.Cancel)
            message.exec_()
            clicked = message.clickedButton()
            if clicked not in (select_button, create_button):
                return None

            creating = clicked == create_button
            dialog = QFileDialog.getSaveFileName if creating else QFileDialog.getOpenFileName
            selected = dialog(
                parent, 'Create configuration' if creating else 'Select configuration',
                '', 'Configuration files (*.json)'
            )
            # Support both iolite's path result and bindings returning (path, filter).
            path = selected[0] if isinstance(selected, (tuple, list)) else selected
            if not path:
                continue
            if creating and not path.lower().endswith('.json'):
                path += '.json'
            path = os.path.abspath(path)
            if creating and os.path.exists(path):
                if QMessageBox.question(
                    parent, 'Replace configuration?',
                    f'{path} already exists. Replace it with the shipped defaults?',
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No
                ) != QMessageBox.Yes:
                    continue
            try:
                if creating:
                    self.write(path, self.defaults)
                config = self.load(path)
            except (OSError, ValueError) as error:
                QMessageBox.warning(parent, 'Configuration error', f'{path}\n\n{error}')
                continue
            self.remember(path, parent)
            return path, config

    def save_updates(self, path, updates, parent=None):
        """Return the saved configuration, or None on cancellation or failure."""
        try:
            overrides = self.read_overrides(path)
            overrides.update(updates)
            config = self.validate(overrides)
            if QMessageBox.question(
                parent, 'Save defaults?',
                'Update the active configuration with the current settings? '
                'This affects everyone using this file.',
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            ) != QMessageBox.Yes:
                return None
            self.write(path, overrides)
        except (OSError, ValueError) as error:
            QMessageBox.warning(parent, 'Could not save defaults', str(error))
            return None
        QMessageBox.information(parent, 'Defaults saved', f'Saved defaults to {path}')
        return config
