"""Ordered registries used to plug new physics into the evolution.

Neutrino processes, thermal-conductivity channels and heat-capacity terms are
all collected in a ``Registry``. Adding or testing a physical ingredient means
registering one function; replacing an existing one means registering under
the same name with ``replace=True``. Nothing else needs to be edited, as long
as the module that registers the function is imported before the run.

Registered functions are evaluated, and their results summed, in registration
order. ``replace=True`` keeps the position of the replaced function, whereas
removing and registering again moves it to the end (this changes only the
round-off of the sums, at the 1e-15 level).
"""


class Registry:
    """Ordered collection of named functions of one kind."""

    def __init__(self, kind):
        self.kind = kind
        self._functions = {}

    def register(self, name, *, replace=False):
        """Decorator: register a function under ``name``.

        An existing name is an error unless ``replace=True``, in which case the
        function takes over the position of the one it replaces.
        """

        def decorator(function):
            if name in self._functions and not replace:
                raise ValueError(
                    f"{self.kind} '{name}' is already registered; use replace=True to override it"
                )
            self._functions[name] = function
            return function

        return decorator

    def unregister(self, name):
        """Remove and return a registered function (for example to switch a process off)."""
        if name not in self._functions:
            raise KeyError(f"No {self.kind} named '{name}'; registered: {self.names()}")
        return self._functions.pop(name)

    def names(self):
        return tuple(self._functions)

    def items(self):
        return tuple(self._functions.items())

    def restore(self, items):
        """Restore a snapshot returned by items(), including its original order."""
        self._functions = dict(items)

    def __contains__(self, name):
        return name in self._functions

    def __len__(self):
        return len(self._functions)
