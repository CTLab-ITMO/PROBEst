# Contributing to PROBESt

***Warning! Tool is under active development***

Thank you for your interest in contributing to PROBESt! We welcome contributions from the community to help improve and expand the functionality of this tool. Whether you're fixing bugs, adding new features, or improving documentation, your contributions are highly appreciated.

This guide provides an overview of how you can contribute to the project. Please take a moment to read through it before making your contribution.

---

## Table of Contents
1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
   - [Setting Up the Development Environment](#setting-up-the-development-environment)
   - [Forking the Repository](#forking-the-repository)
3. [Making Contributions](#making-contributions)
   - [Reporting Bugs](#reporting-bugs)
   - [Suggesting Enhancements](#suggesting-enhancements)
   - [Submitting Pull Requests](#submitting-pull-requests)
4. [Coding Guidelines](#coding-guidelines)
5. [Testing](#testing)
6. [Documentation](#documentation)
7. [Code Review Process](#code-review-process)
8. [Contact](#contact)

---

## Code of Conduct

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md). Please ensure that your contributions align with the values and guidelines outlined in the document.

---

## Getting Started

### Setting Up the Development Environment

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/CTLab-ITMO/PROBEst.git
   cd PROBESt
   ```

2. **Install Dependencies**:
   Ensure you have Python 3.x installed. Then, install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. **Set Up Bash Scripts**:
   Ensure the bash scripts in the `scripts/` directory are executable:
   ```bash
   chmod +x scripts/*.sh
   ```

4. **Run Tests**:
   Verify that the tests pass:
   ```bash
   python -m unittest discover tests/
   ```

### Forking the Repository

1. Fork the [PROBESt repository](https://github.com/CTLab-ITMO/PROBEst) to your GitHub account.
2. Clone your forked repository to your local machine:
   ```bash
   git clone https://github.com/CTLab-ITMO/PROBEst.git
   cd PROBESt
   ```
3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/CTLab-ITMO/PROBEst
   ```

---

## Making Contributions

### Reporting Bugs

If you encounter a bug, please follow these steps:

1. Check the [Issues](https://github.com/CTLab-ITMO/PROBEst/issues) to see if the bug has already been reported.
2. If not, open a new issue with the following details:
   - A clear and descriptive title.
   - Steps to reproduce the bug.
   - Expected vs. actual behavior.
   - Screenshots or error logs (if applicable).

### Suggesting Enhancements

If you have an idea for a new feature or improvement:

1. Check the [Issues](https://github.com/CTLab-ITMO/PROBEst/issues) to see if it has already been suggested.
2. If not, open a new issue with the following details:
   - A clear and descriptive title.
   - A detailed explanation of the enhancement.
   - Use cases or examples (if applicable).

### Submitting Pull Requests

1. Create a new branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. Make your changes and commit them with clear and concise messages:
   ```bash
   git commit -m "Add: Description of your changes"
   ```
3. Push your changes to your forked repository:
   ```bash
   git push origin feature/your-feature-name
   ```
4. Open a Pull Request (PR) against the `main` branch of the original repository. Include:
   - A clear description of the changes.
   - Reference to any related issues.
   - Screenshots or test results (if applicable).

---

## Coding Guidelines

- Follow the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide for Python code.
- Use descriptive variable and function names.
- Add comments to explain complex logic or algorithms.
- Ensure your code is well-documented and includes unit tests.

---

## Testing

- Write unit tests for new features or bug fixes.
- Ensure all tests pass before submitting a PR:
  ```bash
  python -m unittest discover tests/
  ```
- Add integration tests for new scripts or workflows.

---

## Documentation

- Update the documentation (e.g., README, CONTRIBUTING, or inline comments) to reflect your changes.
- If you add a new feature, provide usage examples in the documentation.

---

## Code Review Process

1. Your PR will be reviewed by the maintainers.
2. Address any feedback or requested changes.
3. Once approved, your changes will be merged into the `main` branch.

---

## Contact

If you have any questions or need assistance, feel free to reach out:

- **Email**: dvsmutin@itmo.ru
- **GitHub Issues**: [Open an Issue](https://github.com/CTLab-ITMO/PROBEst/issues)

---

Thank you for contributing to PROBESt! Your efforts help make this tool better for everyone. 🚀