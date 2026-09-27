"""Classed errors: every failure of the video pipeline says what is wrong."""
from __future__ import annotations


class VideoStageError(Exception):
    """Base class for problems the video pipeline detects."""


class InputError(VideoStageError):
    """A pipeline output, audio file or tool this stage needs is missing or malformed."""


class WalkError(VideoStageError):
    """The walk is not a valid path through the joint network."""


class ModelError(VideoStageError):
    """The model could not be called, or did not return valid lyrics."""


class RenderError(VideoStageError):
    """The fly-through could not be rendered or encoded."""
