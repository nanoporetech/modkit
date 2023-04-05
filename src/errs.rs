use std::error::Error;
use std::fmt::Formatter;

#[derive(Debug)]
pub struct InputError(pub String);

impl InputError {
    pub fn new(err: &str) -> Self {
        Self(err.to_owned())
    }
}

impl std::fmt::Display for InputError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl Error for InputError {}

impl From<InputError> for String {
    fn from(err: InputError) -> Self {
        err.0
    }
}

impl From<String> for InputError {
    fn from(s: String) -> Self {
        Self(s)
    }
}

impl From<&str> for InputError {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}

pub enum RunError {
    BadInput(InputError),
    Skipped(String),
    Failed(String),
}

impl RunError {
    pub fn new_input_error<T: Into<String>>(reason: T) -> Self {
        Self::BadInput(InputError(reason.into()))
    }

    pub fn new_skipped<T: Into<String>>(reason: T) -> Self {
        Self::Skipped(reason.into())
    }

    pub fn new_failed<T: Into<String>>(reason: T) -> Self {
        Self::Failed(reason.into())
    }
}

impl From<InputError> for RunError {
    fn from(err: InputError) -> Self {
        Self::BadInput(err)
    }
}

impl std::fmt::Display for RunError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let reason = match self {
            RunError::BadInput(InputError(reason)) => {
                format!("Bad Input: {}", reason)
            }
            RunError::Skipped(reason) => format!("Skipped: {}", reason),
            RunError::Failed(reason) => format!("Failed: {}", reason),
        };
        write!(f, "{}", reason)
    }
}

impl std::fmt::Debug for RunError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl Error for RunError {}
