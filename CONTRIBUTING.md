# Contributing to BZea_Genotype

Thank you for your interest in contributing to this ANGSD genotyping tutorial for maize! We welcome contributions from the community.

## How to Contribute

### Reporting Issues

If you find bugs, errors in documentation, or have suggestions:

1. Check if the issue already exists in [GitHub Issues](../../issues)
2. If not, create a new issue with:
   - Clear, descriptive title
   - Detailed description of the problem
   - Steps to reproduce (if applicable)
   - Your environment (OS, ANGSD version, etc.)
   - Expected vs actual behavior

### Suggesting Enhancements

For new features or improvements:

1. Open an issue describing your suggestion
2. Explain the use case and why it would be valuable
3. Discuss the proposed implementation

### Contributing Code or Documentation

1. **Fork the repository**
   ```bash
   # Click "Fork" on GitHub, then clone your fork
   git clone https://github.com/YOUR_USERNAME/BZea_Genotype.git
   cd BZea_Genotype
   ```

2. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bug-fix
   ```

3. **Make your changes**
   - Follow existing code style and documentation format
   - Test your changes thoroughly
   - Add comments where necessary
   - Update documentation if needed

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "Add: brief description of your changes"
   ```

5. **Push and create a Pull Request**
   ```bash
   git push origin feature/your-feature-name
   ```
   Then open a Pull Request on GitHub

## Contribution Guidelines

### Code Style

**Bash scripts:**
- Use descriptive variable names (UPPERCASE for constants)
- Add comments explaining complex logic
- Include usage examples in header
- Check scripts with `shellcheck` if possible

**R scripts:**
- Follow tidyverse style guide
- Use meaningful variable names
- Comment complex operations
- Test with example data

### Documentation

- Use clear, concise language
- Include code examples
- Provide expected outputs
- Keep tutorials beginner-friendly
- Use proper markdown formatting

### Adding New Scripts

When contributing new analysis scripts:

1. Add comprehensive header with:
   - Purpose of the script
   - Usage instructions
   - Required inputs
   - Expected outputs
   - Author information

2. Include error checking:
   ```bash
   if [ ! -f "${REQUIRED_FILE}" ]; then
       echo "ERROR: Required file not found"
       exit 1
   fi
   ```

3. Make scripts configurable with command-line arguments

4. Test with example data

5. Add documentation to `scripts/README.md`

### Adding New Documentation

When adding or updating documentation:

1. Use clear section headings
2. Include practical examples
3. Add troubleshooting sections
4. Link to related documentation
5. Test all commands provided

### Example Contributions We'd Love

- **Additional analysis scripts**
  - Admixture analysis pipeline
  - Selection scan methods
  - Visualization improvements
  
- **Documentation improvements**
  - More detailed examples
  - Common error solutions
  - Best practices for specific use cases
  
- **Tutorial enhancements**
  - Video tutorials
  - Jupyter notebooks
  - Additional worked examples
  
- **Data resources**
  - Links to public datasets
  - Example analysis outputs
  - Benchmark results

## Testing

Before submitting:

1. **Test your scripts** with example data
2. **Check documentation** for accuracy
3. **Verify links** work correctly
4. **Run spell check** on documentation
5. **Ensure code follows** style guidelines

## Code of Conduct

### Our Standards

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on what's best for the community
- Show empathy towards others

### Unacceptable Behavior

- Harassment or discriminatory language
- Personal attacks
- Trolling or insulting comments
- Publishing others' private information
- Other unethical or unprofessional conduct

## Questions?

- Open an issue for questions about contributing
- Reach out to maintainers for guidance
- Check existing issues and PRs for similar contributions

## Recognition

Contributors will be acknowledged in:
- Repository contributors list
- Release notes (for significant contributions)
- Documentation (for major additions)

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (see LICENSE file).

## Getting Started Checklist

- [ ] Read this contributing guide
- [ ] Check existing issues and PRs
- [ ] Fork the repository
- [ ] Create a feature branch
- [ ] Make your changes
- [ ] Test thoroughly
- [ ] Update documentation
- [ ] Commit with clear messages
- [ ] Submit Pull Request
- [ ] Respond to review feedback

Thank you for helping make this resource better for the maize genetics community! 🌽
