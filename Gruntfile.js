module.exports = function(grunt) {
  grunt.initConfig({
    pkg: grunt.file.readJSON('package.json'),
    jshint: {
      all: ['<Grunt>  </Grunt>file.js',
            'src/**/*.js']
    },
    browserify: {
      production: {
        src: 'src/program.js',
        dest: 'build/program.js',
        options: {
          externalize: ['src/program.js']
        }
      },
      debug: {
        src: 'src/program.js',
        dest: 'build/program.debug.js',
        options: {
          externalize: ['src/program.js'],
          debug: true
        }
      }
    },
    uglify: {
      build: {
        src: 'build/program.js',
        dest: 'build/program.min.js'
      }
    },
    'string-replace': {
      kit: {
        files: {
          'build/program.min.js' : 'build/program.min.js',
          'build/program.js' : 'build/program.js',
        },
        options: {
          replacements: [{
            pattern: /\brequire\b/ig,
            replacement: 'req'
          }]
        }
      }
    },
    watch: {
      scripts: {
        files: ['src/**/*.js'],
        tasks: ['jshint',
                'browserify:production',
                'browserify:debug',
                 'uglify',
                 'string-replace'],
        options: {
          nospawn: true
        }
      }
    }
  });

  grunt.loadNpmTasks('grunt-contrib-uglify');
  grunt.loadNpmTasks('grunt-contrib-jshint');
  grunt.loadNpmTasks('grunt-browserify');
  grunt.loadNpmTasks('grunt-contrib-watch');
  grunt.loadNpmTasks('grunt-string-replace');


  grunt.registerTask('default', ['jshint',
                                 'browserify:production',
                                 'uglify',
                                 'string-replace']);

  grunt.registerTask('debug', ['jshint',
                               'browserify:debug']);

};